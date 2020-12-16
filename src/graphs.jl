# Each pixel is given a unique ID, missings are automatically considered no data
# and elements with value equal to no_data_val are also ignored (not considered valid nodes)
function construct_nodemap(weights::Matrix{T} where T <: Real;
                           no_data_val = nothing)
    dims = size(weights)

    # Make an resistance of unique node identifiers
    nodemap = zeros(Int64, dims)
    is_node = coalesce.(weights .!= no_data_val, false)
    nodemap[is_node] = 1:sum(is_node)

    nodemap
end

function construct_nodemap(weights::GeoData.GeoArray)
    # Handle dimensions, get them so the order matches the input
    lat_lon_dims = get_lat_lon_dims(weights)

    # Make an resistance of unique node identifiers
    nodemap = zeros(Int64, size(weights.data[:, :, 1]))
    is_node = weights.data[:, :, 1] .!= weights.missingval
    nodemap[is_node] = 1:sum(is_node)

    nodemap = GeoData.GeoArray(nodemap, dims = lat_lon_dims)

    nodemap
end

# Add GeoArray method
function construct_graph(weights::Matrix{T} where T <: Number,
                         nodemap::Matrix{Int};
                         no_data_val = nothing,
                         weights_layer_is_conductance::Bool = false,
                         connect_four_neighbors_only::Bool = false,
                         connect_using_avg_weights::Bool = true)
    # Which averaging function to use
    card_avg = connect_using_avg_weights ? res_cardinal_avg : cond_cardinal_avg
    diag_avg = connect_using_avg_weights ? res_diagonal_avg : cond_diagonal_avg
    dims = size(weights)
    not_no_data = weights .!= no_data_val

    if sum((weights .<= 0) .& not_no_data) != 0
        @error "weights contains 0 or negative values aside from the provided no_data_val, which is not supported" && return
    end

    if weights_layer_is_conductance
        weights[not_no_data] = 1 ./ weights[not_no_data] # TODO I think this is mutating the input somehow... double check
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
                    res = card_avg(weights[row, column],
                                   weights[row, column + 1])
                    push!(sources, nodemap[row, column])
                    push!(destinations, nodemap[row, column + 1])
                    push!(node_weights, res)
                end
                # South
                if row != dims[1] && nodemap[row + 1, column] != 0
                    res = card_avg(weights[row, column],
                                   weights[row + 1, column])
                    push!(sources, nodemap[row, column])
                    push!(destinations, nodemap[row + 1, column])
                    push!(node_weights, res)
                end

                ## Add diagonal neighbors if needed
                if !connect_four_neighbors_only
                    # Northeast
                    if column != dims[2] && row != 1 && nodemap[row - 1, column + 1] != 0
                        res = diag_avg(weights[row, column],
                                       weights[row - 1, column + 1])
                        push!(sources, nodemap[row, column])
                        push!(destinations, nodemap[row - 1, column + 1])
                        push!(node_weights, res)
                    end
                    # Southeast
                    if row != dims[1] && column != dims[2] && nodemap[row + 1, column + 1] != 0
                        res = diag_avg(weights[row, column],
                                       weights[row + 1, column + 1])
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

function construct_graph(weights::GeoData.GeoArray,
                         nodemap::GeoData.GeoArray;
                         weights_layer_is_conductance::Bool = false,
                         connect_four_neighbors_only::Bool = false,
                         connect_using_avg_weights::Bool = true)
    construct_graph(weights.data[:, :, 1],
                    nodemap.data[:, :, 1],
                    no_data_val = weights.missingval,
                    weights_layer_is_conductance = weights_layer_is_conductance,
                    connect_four_neighbors_only = connect_four_neighbors_only,
                    connect_using_avg_weights = connect_using_avg_weights)
end
