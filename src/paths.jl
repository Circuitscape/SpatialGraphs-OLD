# TODO add a function that takes output from this and converts it to a GeoArray or Matrix
function compute_cd_graph(g::AbstractSimpleWeightedGraph,
                          node_ids; # will report lowest CD to get to any of these Number or Vector
                          dist_fun::Function = dijkstra_shortest_paths)
    cd_graph = dist_fun(g, node_ids)

    return cd_graph
end

function sample_lcp_node_pairs(sample_weights::Matrix{T} where T <: Number, # Matrix of weights
                               nodemap::Matrix{Int}, # Matrix of node_ids to sample
                               n_pairs::Int)
    # Set weights (habitat) to 0 for non-nodes
    sample_weights[nodemap .== 0] .= 0

    # Mask weights objects from habitat layer for use in sampling
    sample_weights_weights = StatsBase.weights(sample_weights)

    # Initialize object in which to store pairs of nodemap
    samples = Vector{Vector{Int}}()

    # Sample n_pairs * 2 points
    the_sample = sample(nodemap, sample_weights_weights, n_pairs * 2, replace = true)

    for i in 1:(n_pairs)
        push!(samples, the_sample[(2 * i - 1):(2 * i)])
    end

    return samples
end

function sample_lcp_node_pairs(sample_weights::GeoData.GeoArray,
                               nodemap::GeoData.GeoArray,
                               n_pairs::Int)
    weights = deepcopy(sample_weights.data)
    weights[weights .== sample_weights.missingval] .= 0
    samples = sample_lcp_node_pairs(sample_weights.data[:, :, 1],
                                    nodemap.data[:, :, 1],
                                    n_pairs)

    return samples
end

function least_cost_path(g::AbstractGraph, start, destination::Int)
    cost = compute_cd_graph(g, start)
    return enumerate_paths(cost, destination)
end

# Returns a nodemap to map node_ids to points in space
# and a vector of vectors (each vector is a path)
# TODO update this to use a Matrix or GeoData.GeoArray
function random_lcps(resistance::Array{Union{Missing, T}, 2} where T <: Number,
                     weights::Array{Union{Missing, T}, 2} where T <: Number,
                     n_paths::Int;
                     source_percentile_threshold::Number = 0,
                     resistance_layer_is_conductance::Bool = false,
                     connect_four_neighbors_only::Bool = false,
                     connect_using_avg_resistances::Bool = false,
                     parallel::Bool = true)
    if (source_percentile_threshold > 0)
        weights = threshold_array!(weights, source_percentile_threshold)
    end

    nodemap = construct_nodemap(resistance)
    graph = construct_graph(resistance,
                            nodemap,
                            resistance_layer_is_conductance = resistance_layer_is_conductance,
                            connect_four_neighbors_only = connect_four_neighbors_only,
                            connect_using_avg_resistances = connect_using_avg_resistances)

    node_pairs = sample_lcp_node_pairs(weights,
                                       nodemap,
                                       n_pairs)
    #### Identify the least cost paths
    # Initialize the paths
    paths = Vector{Vector{Int}}(undef, n_pairs)

    if parallel
        @threads for i in 1:n_pairs
           lcp = least_cost_path(graph, node_pairs[i][1], node_pairs[i][2])
           paths[i] = lcp
        end
    else
        for i in 1:n_pairs
           lcp = least_cost_path(graph, node_pairs[i][1], node_pairs[i][2])
           paths[i] = lcp
        end
    end

    return nodemap, paths
end

# TODO add method that accepts GeoData.GeoArray nodemap
function path_to_array(path::Vector{Int},
                       nodemap::Matrix{Int})
    path_array = zeros(eltype(nodemap), size(nodemap))
    path_array[in.(nodemap, [path])] .= 1

    return path_array
end

# Convert a path to a vector of it's coordinates, if the nodemap
# TODO add method that accepts GeoData.GeoArray nodemap
function path_to_points(path::Vector{Int},
                        nodemap::Matrix{Int},
                        geotransform::Vector{N} where N <: Number;
                        parallel::Bool = true)
    cart_coords = Vector{Tuple{Int64, Int64}}(undef, length(path))

    if parallel
        @threads for i in 1:length(path)
            cart_coord = findall(nodemap .== path[i])
            cart_coords[i] = convert.(Tuple, cart_coord)[1]
        end
    else
        for i in 1:length(path)
            cart_coord = findall(nodemap .== path[i])
            cart_coords[i] = convert.(Tuple, cart_coord)[1]
        end
    end

    upper_left_corner_coords = (geotransform[1], geotransform[4])
    multiplier = (geotransform[2], geotransform[6])

    # Align points with pixel center instead of upper left corners
    upper_left_center_coords = upper_left_corner_coords .+ 0.5 .* multiplier

    # Map cartesian coordinate to geographic coordinate
    geo_coords = map(x -> x .* multiplier .+ upper_left_center_coords, cart_coords)

    geo_coords
end
