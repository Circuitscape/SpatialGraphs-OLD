using Test, SpatialGraphs, LightGraphs, GeoData, ArchGDAL

url_base = "https://raw.githubusercontent.com/Circuitscape/datasets/main/"
# Download the NLCD tile used to create the resistance surface and load it
download(string(url_base, "data/nlcd_2016_frederick_md.tif"),
         "nlcd_2016_frederick_md.tif")

garray = GeoArray(GDALarray("nlcd_2016_frederick_md.tif", missingval = -9))

@testset "graph construction" begin
    no_data_val = -9999
    weights = [1 1 1;
               1 3 1;
               1 NaN no_data_val]

    nodemap = construct_nodemap(weights, no_data_val = no_data_val)

    # NoData entries in weights are 0 in nodemap
    @test (nodemap .== 0) == ((weights .== no_data_val) .| isnan.(weights))

    # Values in nodemap are as expected
    @test sort(collect(nodemap[nodemap .!= 0])) == collect(1:sum((weights .!= no_data_val) .& (!).(isnan.(weights))))

    ## Check that values in the resistance graph are as expected
    graph = construct_graph(weights, nodemap, no_data_val = no_data_val)
    the_edges = collect(edges(graph))

    # Test that the edges are correct and have proper weights
    for i in 1:length(the_edges)
        source_i = src(the_edges[i])
        dest_i = dst(the_edges[i])
        weight_i = the_edges[i].weight

        source_coords = findall(nodemap .== source_i)[1]
        dest_coords = findall(nodemap .== dest_i)[1]

        # Test that source row is within 1 step of dest row
        row_diff = abs(source_coords[1] - dest_coords[1])
        @test row_diff <= 1

        # Test that source column is within 1 step of dest row
        col_diff = abs(source_coords[2] - dest_coords[2])
        @test col_diff <= 1

        # Test that the weight is what it should be (assumes connect_using_avg_resistance = true in graph construction)
        if (row_diff == 1 && col_diff == 1) # get diagonal average
            @test weight_i == SpatialGraphs.res_diagonal_avg(weights[source_coords], weights[dest_coords])
        else
            @test weight_i == SpatialGraphs.res_cardinal_avg(weights[source_coords], weights[dest_coords])
        end
    end

    ## GeoArray compatibility
    old_data = deepcopy(garray.data[:, :, 1])
    nodemap = construct_nodemap(garray)

    @test size(nodemap) == size(garray[Band(Between(1, 1))])
    g = construct_graph(garray, nodemap,
                        cost_layer_is_conductance = true,
                        connect_four_neighbors_only = true,
                        connect_using_avg_resistance = false)

    # Make sure none of the functions changed input
    @test old_data == float.(garray.data[:, :, 1])
end

@testset "internals" begin
    # get_lat_lon_dims, test that they are returned in the correct order
    ll_dims = SpatialGraphs.get_lat_lon_dims(garray)
    @test ll_dims[1] == dims(garray)[1]
    @test ll_dims[2] == dims(garray)[2]
end

# LightGraphs and SimpleWeightedGraphs presumably already test for correctness of
# methods, so above tests to confirm graphs are constructed properly are most
# important
@testset "paths" begin
    no_data_val = -9999
    weight = [1 1 1;
              1 3 1;
              1 1 1]

    nodemap = construct_nodemap(weight, no_data_val = no_data_val)
    g = construct_graph(weight, nodemap, no_data_val = no_data_val)

    cdist_array = cost_distance(g, nodemap, 5)
    @test size(cdist_array) == size(nodemap)
    @test all(cdist_array[nodemap[nodemap .∈ [collect(2:2:8)]]] .== 2.0)
    @test all(cdist_array[nodemap[nodemap .∈ [[1, 3, 7]]]] .==
        SpatialGraphs.res_diagonal_avg(3, 1))

    lcps, nodemap = random_lcps(weight, fill(1., size(weight)), 50, parallel = false)

    path = least_cost_path(g, 2, 8)
    path_coords = path_to_cartesian_coords(path, nodemap)
    path_array = path_to_array(path, nodemap)

    # GeoArray stuff
    nodemap = construct_nodemap(garray)
    g = construct_graph(garray, nodemap,
                        cost_layer_is_conductance = true,
                        connect_four_neighbors_only = false,
                        connect_using_avg_resistance = false)
    path = least_cost_path(g, 2, 8)
    path_geoarray = path_to_geoarray(path, nodemap)
    path_linestring = path_to_linestring(path, nodemap)

    cd_geoarray = cost_distance(g, nodemap, [2, 3])

    # random lcps, and convert to geometries
    lcps, nodemap = random_lcps(garray, garray, 10)
    linestrings = paths_to_linestrings(lcps, nodemap)
end

rm("nlcd_2016_frederick_md.tif")

