using Test, SpatialGraphs, LightGraphs

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
