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
               1 1 no_data_val]

    nodemap = construct_nodemap(weights, no_data_val = no_data_val)

    # NoData entries in weights are 0 in nodemap
    @test (nodemap .== 0) == (weights .== no_data_val)

    # Values in nodemap are as expected
    @test sort(collect(nodemap[nodemap .!= 0])) == collect(1:sum(weights .!= no_data_val))

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
end

@testset "GeoArray compatibility" begin
    old_data = deepcopy(garray.data[:, :, 1])

    nodemap = construct_nodemap(garray)

    @test size(nodemap) == size(garray.data[:, :, 1])
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

rm("nlcd_2016_frederick_md.tif")
