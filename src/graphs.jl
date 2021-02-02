# Each pixel is given a unique ID, missings are automatically considered no data
# and elements with value equal to no_data_val are also ignored (not considered valid nodes)
"""
    construct_nodemap(A::Matrix{T} where T <: Real;
                      no_data_val = nothing)
    construct_nodemap(A::GeoData.GeoArray)

Construct a nodemap from `A` and with dims equal to `size(A)` that can be used
to provide spatial references for each vertex in a graph representation of
`A`. Each element in the nodemap is given a unique integer ID. Elements in the
nodemap that correspond to NoData in `A` (`A.missingval` if `A` is a GeoArray,
or `no_data_val` if `A` is a matrix) are given a value of 0.
"""
function construct_nodemap(A::Matrix{T} where T<: Real;
                           no_data_val = nothing)
    dims = size(A)

    # Make an resistance of unique node identifiers
    nodemap = zeros(Int64, dims)
    is_node = (A .!= no_data_val) .&
        ((!).(isnan.(A)))
    nodemap[is_node] = 1:sum(is_node)

    nodemap
end

function construct_nodemap(A::GeoData.GeoArray)
    # Handle dimensions, get them so the order matches the input
    lat_lon_dims = get_lat_lon_dims(A)
    A_band1 = A[Band(Between(1, 1))]

    # Make an array of unique node identifiers
    nodemap = zeros(Int64, size(A_band1.data))
    is_node = (A_band1.data .!= A.missingval) .&
        ((!).(isnan.(A_band1.data)))

    nodemap[is_node] = 1:sum(is_node)

    nodemap = GeoData.GeoArray(nodemap, dims = lat_lon_dims)

    nodemap
end

"""
    construct_weighted_graph(cost_surface::Matrix{T} where T <: Real,
                             nodemap::Matrix{Int};
                             no_data_val = nothing,
                             cost_layer_is_conductance = false,
                             connect_four_neighbors_only = false,
                             connect_using_avg_cost = true)

    construct_weighted_graph(cost_surface::GeoData.GeoArray,
                             nodemap::GeoData.GeoArray;
                             cost_layer_is_conductance = false,
                             connect_four_neighbors_only = false,
                             connect_using_avg_cost = true)

Construct a impleWeightedGraph from a `cost_surface` representing traversal costs`
and a `nodemap` representing the locations in space (either geographic or
cartesian) of the vertices to be added to the graph. The resulting graph will
contain edges connecting neighboring elements in `nodemap` with weights
calculated using the corresponding elements in `cost_surface` for each element
in the nodemap. For cardinal neighbors with cartesian coordinates [i, j] and
[i+1, j], the corresponding edge connecting node `nodemap[i, j]` to node
`nodemap[i + 1, j]`. It's weight (traversal cost) equals
`(cost_surface[i, j] + cost_surface[i+1, j]) / 2`. For diagonal neighbors,
[i, j] and [i+1, j+1], the edge weight equals
`(2/√2)*(cost_surface[i, j] + cost_surface[i+1, j+1]) / 2`. The average of the
two costs is multiplied by 2/√2 to account for the increased distance between
diagonal neighbors.

## Keyword Arguments
- `no_data_val`: Real. The value corresponding to "no data" in `cost_surface`.
Elements equal to `no_data_val` correspond to pixels that cannot be traversed.
Defaults to `nothing`, or, if cost_surface is a GeoData.GeoArray,
cost_surface.missingval is used.

- `cost_layer_is_conductace`: Boolean. Does the `cost_surface` represent
permeability/conductance instead of cost/resistance? If `true`, cost is
calculated as `1 / cost_surface`. Defaults to `false`.

- `connect_four_neighbors_only`: Boolean. Connect only cardinal neighbors in
`nodemap`? Defaults to `false`.

- `connect_using_avg_cost`: `Boolean`. This is intended to offer methods that
complement those used in Circuitscape.jl and Omniscape.jl. In this context,
cost is in units of electrical resistance. If `false`, the cost between two
nodes with resistances R1 and R2 is calculated by converting resistance to
conductances, taking the average, then taking the inverse of the result to
convert back to resistance: `1 / ((1/R1 + 1/R2) / 2)`.
`connect_using_avg_cost = false` correspondes to the default settings in
Circuitscape. Defaults to `true'.
"""
function construct_weighted_graph(cost_surface::Matrix{T} where T <: Real,
                                  nodemap::Matrix{Int};
                                  no_data_val = nothing,
                                  cost_layer_is_conductance::Bool = false,
                                  connect_four_neighbors_only::Bool = false,
                                  connect_using_avg_cost::Bool = true)
    cost = float.(deepcopy(cost_surface))

    # Which averaging function to use
    card_avg = connect_using_avg_cost ? res_cardinal_avg : cond_cardinal_avg
    diag_avg = connect_using_avg_cost ? res_diagonal_avg : cond_diagonal_avg
    dims = size(cost)
    not_no_data = cost .!= no_data_val

    if sum((cost .<= 0) .& not_no_data) != 0
        @error("cost surface contains 0 or negative values (aside from the provided no_data_val), which is not supported")
    end

    if cost_layer_is_conductance
        cost[not_no_data] = 1 ./ cost[not_no_data]
    end

    # Construct graph
    #res_graph = SimpleWeightedGraph(sum(is_node))
    sources = Vector{Int64}()
    destinations = Vector{Int64}()
    node_weights = Vector{Float64}()

    # Add the edges
    # only need to do neighbors down or to the right because
    # the graph is undirected and edge additions will be redundant
    for row in 1:dims[1]
        for column in 1:dims[2]
            if nodemap[row, column] == 0 || cost[row, column] == no_data_val
                continue
            else
                ## Add cardinal neighbors
                # East
                if column != dims[2] && nodemap[row, column + 1] != 0 && cost[row, column + 1] != no_data_val
                    res = card_avg(cost[row, column],
                                   cost[row, column + 1])
                    push!(sources, nodemap[row, column])
                    push!(destinations, nodemap[row, column + 1])
                    push!(node_weights, res)
                end
                # South
                if row != dims[1] && nodemap[row + 1, column] != 0  && cost[row + 1, column] != no_data_val
                    res = card_avg(cost[row, column],
                                   cost[row + 1, column])
                    push!(sources, nodemap[row, column])
                    push!(destinations, nodemap[row + 1, column])
                    push!(node_weights, res)
                end

                ## Add diagonal neighbors if needed
                if !connect_four_neighbors_only
                    # Northeast
                    if column != dims[2] && row != 1 && nodemap[row - 1, column + 1] != 0 && cost[row - 1, column + 1] != no_data_val
                        res = diag_avg(cost[row, column],
                                       cost[row - 1, column + 1])
                        push!(sources, nodemap[row, column])
                        push!(destinations, nodemap[row - 1, column + 1])
                        push!(node_weights, res)
                    end
                    # Southeast
                    if row != dims[1] && column != dims[2] && nodemap[row + 1, column + 1] != 0  && cost[row + 1, column + 1] != no_data_val
                        res = diag_avg(cost[row, column],
                                       cost[row + 1, column + 1])
                        push!(sources, nodemap[row, column])
                        push!(destinations, nodemap[row + 1, column + 1])
                        push!(node_weights, res)
                    end
                end
            end
        end
    end

    # push!'ing to vectors then combining is way faster than add_edge!
    # Sometimes, depending on the nodemap, there may be multiple edges with
    # different weights connecting the same two vertices (only when vertices
    # span multiple "pixels"). Default to using the min of all duplicates via
    # combine = min
    resistance_graph = SimpleWeightedGraph(sources, destinations, node_weights, combine = min)

    return resistance_graph
end

function construct_weighted_graph(cost_surface::GeoData.GeoArray,
                                  nodemap::GeoData.GeoArray;
                                  cost_layer_is_conductance::Bool = false,
                                  connect_four_neighbors_only::Bool = false,
                                  connect_using_avg_cost::Bool = true)
    construct_weighted_graph(cost_surface.data[:, :, 1],
                             nodemap.data[:, :, 1],
                             no_data_val = cost_surface.missingval,
                             cost_layer_is_conductance = cost_layer_is_conductance,
                             connect_four_neighbors_only = connect_four_neighbors_only,
                             connect_using_avg_cost = connect_using_avg_cost)
end

function get_cartesian_indices(node_array::Matrix{Int})
    # Get Cartesian Index of all nodes once, this method ensures that the index
    # in all_coords (created below) matches the node value, e.g.
    # all_coords[10] is the cartesian index for the node with value 10
    ###########################################
    # THIS WILL ONLY WORK WITH NODEMAPS THAT ##
    # WERE CREATED WITH construct_nodemap()  ##
    ###########################################

    # Throw an error if node is formatted the right way
    nodes_bool = node_array .!= 0
    if node_array[nodes_bool] != collect(1:sum(nodes_bool))
        @error("This function does not support nodemaps that were not created with construct_nodemap")
    end


    all_coords = Vector{CartesianIndex{2}}(undef, maximum(node_array))
    size_dims = size(node_array)
    idx = 1
    # Need to do cols then rows to match how nodemap is created
    # this is WAY faster than a bunch of findall's
    for col in 1:size_dims[2]
        for row in 1:size_dims[1]
            if node_array[row, col] == 0
                continue
            end
            all_coords[idx] = CartesianIndex(row, col)
            idx+=1
        end
    end

    return all_coords
end

function get_cartesian_indices(nodemap::GeoArray)
    get_cartesian_indices(nodemap.data[:, :, 1])
end
