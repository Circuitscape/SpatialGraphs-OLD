# Each pixel is given a unique ID, missings are automatically considered no data
# and elements with value equal to no_data_val are also ignored (not considered valid nodes)
function construct_nodemap(cost_surface::Matrix{T} where T <: Real;
                           no_data_val = nothing)
    dims = size(cost_surface)

    # Make an resistance of unique node identifiers
    nodemap = zeros(Int64, dims)
    is_node = (cost_surface .!= no_data_val) .&
        ((!).(isnan.(cost_surface)))
    nodemap[is_node] = 1:sum(is_node)

    nodemap
end

function construct_nodemap(cost_surface::GeoData.GeoArray)
    # Handle dimensions, get them so the order matches the input
    lat_lon_dims = get_lat_lon_dims(cost_surface)
    cost_surface_band1 = cost_surface[Band(Between(1, 1))]

    # Make an array of unique node identifiers
    nodemap = zeros(Int64, size(cost_surface_band1.data))
    is_node = (cost_surface_band1.data .!= cost_surface.missingval) .&
        ((!).(isnan.(cost_surface_band1.data)))

    nodemap[is_node] = 1:sum(is_node)

    nodemap = GeoData.GeoArray(nodemap, dims = lat_lon_dims)

    nodemap
end

# Add GeoArray method
function construct_graph(cost_surface::Matrix{T} where T <: Real,
                         nodemap::Matrix{Int};
                         no_data_val = nothing,
                         cost_layer_is_conductance::Bool = false,
                         connect_four_neighbors_only::Bool = false,
                         connect_using_avg_resistance::Bool = true)
    cost = float.(deepcopy(cost_surface))

    # Which averaging function to use
    card_avg = connect_using_avg_resistance ? res_cardinal_avg : cond_cardinal_avg
    diag_avg = connect_using_avg_resistance ? res_diagonal_avg : cond_diagonal_avg
    dims = size(cost)
    not_no_data = cost .!= no_data_val

    if sum((cost .<= 0) .& not_no_data) != 0
        @error "cost surface contains 0 or negative values (aside from the provided no_data_val), which is not supported" && return
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
            if nodemap[row, column] == 0
                continue
            else
                ## Add cardinal neighbors
                # East
                if column != dims[2] && nodemap[row, column + 1] != 0
                    res = card_avg(cost[row, column],
                                   cost[row, column + 1])
                    push!(sources, nodemap[row, column])
                    push!(destinations, nodemap[row, column + 1])
                    push!(node_weights, res)
                end
                # South
                if row != dims[1] && nodemap[row + 1, column] != 0
                    res = card_avg(cost[row, column],
                                   cost[row + 1, column])
                    push!(sources, nodemap[row, column])
                    push!(destinations, nodemap[row + 1, column])
                    push!(node_weights, res)
                end

                ## Add diagonal neighbors if needed
                if !connect_four_neighbors_only
                    # Northeast
                    if column != dims[2] && row != 1 && nodemap[row - 1, column + 1] != 0
                        res = diag_avg(cost[row, column],
                                       cost[row - 1, column + 1])
                        push!(sources, nodemap[row, column])
                        push!(destinations, nodemap[row - 1, column + 1])
                        push!(node_weights, res)
                    end
                    # Southeast
                    if row != dims[1] && column != dims[2] && nodemap[row + 1, column + 1] != 0
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
    resistance_graph = SimpleWeightedGraph(sources, destinations, node_weights)

    return resistance_graph
end

function construct_graph(cost_surface::GeoData.GeoArray,
                         nodemap::GeoData.GeoArray;
                         cost_layer_is_conductance::Bool = false,
                         connect_four_neighbors_only::Bool = false,
                         connect_using_avg_resistance::Bool = true)
    construct_graph(cost_surface.data[:, :, 1],
                    nodemap.data[:, :, 1],
                    no_data_val = cost_surface.missingval,
                    cost_layer_is_conductance = cost_layer_is_conductance,
                    connect_four_neighbors_only = connect_four_neighbors_only,
                    connect_using_avg_resistance = connect_using_avg_resistance)
end

function get_cartesian_indices(node_array::Matrix{Int})
    # Get Cartesian Index of all nodes once, this method ensures that the index
    # in all_coords (created below) matches the node value, e.g.
    # all_coords[10] is the cartesian index for the node with value 10
    ###########################################
    # THIS WILL ONLY WORK WITH NODEMAPS THAT ##
    # WERE CREATED WITH construct_nodemap()  ##
    ###########################################
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

