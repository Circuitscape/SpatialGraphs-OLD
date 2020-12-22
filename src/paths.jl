# TODO add a function that takes output from this and converts it to a GeoArray or Matrix
function cost_distance(g::AbstractSimpleWeightedGraph,
                       node_ids; # will report lowest CD to get to any of these Number or Vector
                       dist_fun::Function = dijkstra_shortest_paths)
    cd_graph = dist_fun(g, node_ids)

    return cd_graph
end

function cd_to_array(cdist::LightGraphs.AbstractPathState,
                     nodemap::Matrix{Int};
                     cost_threshold::Real = Inf)
    cd_array = fill(0., size(nodemap))
    cd_array .= cdist.dists[nodemap]

    cd_array[cd_array .> cost_threshold] .= -1

    cd_array
end

function cd_to_geoarray(cdist::LightGraphs.AbstractPathState,
                        nodemap::GeoData.GeoArray;
                        cost_threshold::Real = Inf)
    cd_array = cd_to_array(cdist,
                           nodemap.data[:, :, 1],
                           cost_threshold = cost_threshold)
    lat_lon_dims = get_lat_lon_dims(nodemap)

    cd_geoarray = GeoData.GeoArray(cd_array,
                                   dims = lat_lon_dims,
                                   missingval = -1)
    cd_geoarray
end

function sample_lcp_node_pairs(sample_weights::Matrix{T} where T <: Real, # Matrix of weights
                               nodemap::Matrix{Int}, # Matrix of node_ids to sample
                               n_pairs::Int)
    weights = deepcopy(sample_weights)
    # Set weights (habitat) to 0 for non-nodes
    weights[nodemap .== 0] .= 0

    # Mask weights objects from habitat layer for use in sampling
    weights_weights = StatsBase.weights(weights)

    # Initialize object in which to store pairs of nodemap
    samples = Vector{Vector{Int}}()

    # Sample n_pairs * 2 points
    the_sample = sample(nodemap, weights_weights, n_pairs * 2, replace = true)

    for i in 1:(n_pairs)
        push!(samples, the_sample[(2 * i - 1):(2 * i)])
    end

    return samples
end

function sample_lcp_node_pairs(sample_weights::GeoData.GeoArray,
                               nodemap::GeoData.GeoArray,
                               n_pairs::Int)
    weights = deepcopy(sample_weights.data[:, :, 1])
    weights[weights .== sample_weights.missingval] .= 0
    samples = sample_lcp_node_pairs(weights,
                                    nodemap.data[:, :, 1],
                                    n_pairs)

    return samples
end

function least_cost_path(g::AbstractGraph,
                         start::Union{Int, Vector{T} where T <: Int},
                         destination::Int)
    cost = compute_cd_graph(g, start)
    return enumerate_paths(cost, destination)
end

# Returns a nodemap to map node_ids to points in space
# and a vector of vectors (each vector is a path)
function random_lcps(weights::Matrix{T} where T <: Real,
                     sample_weights::Matrix{T} where T <: Real,
                     n_paths::Int;
                     no_data_val = nothing,
                     resistance_layer_is_conductance::Bool = false,
                     connect_four_neighbors_only::Bool = false,
                     connect_using_avg_resistances::Bool = false,
                     parallel::Bool = true)
    nodemap = construct_nodemap(weights)
    graph = construct_graph(weights,
                            nodemap,
                            no_data_val = no_data_val,
                            resistance_layer_is_conductance = resistance_layer_is_conductance,
                            connect_four_neighbors_only = connect_four_neighbors_only,
                            connect_using_avg_resistances = connect_using_avg_resistances)

    node_pairs = sample_lcp_node_pairs(sample_weights,
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

function random_lcps(weights::GeoData.GeoArray,
                     sample_weights::GeoData.GeoArray,
                     n_paths::Int;
                     resistance_layer_is_conductance::Bool = false,
                     connect_four_neighbors_only::Bool = false,
                     connect_using_avg_resistances::Bool = false,
                     parallel::Bool = true)
    nodemap, paths = random_lcps(weights.data[:, :, 1],
                                 sample_weights.data[:, :, 1],
                                 n_paths;
                                 no_data_val = weights.missingval,
                                 resistance_layer_is_conductance = resistance_layer_is_conductance,
                                 connect_four_neighbors_only = connect_four_neighbors_only,
                                 connect_using_avg_resistances = connect_using_avg_resistances,
                                 parallel = parallel)
    # Convert nodemap to GeoArray
    lat_lon_dims = get_lat_lon_dims(weights)
    nodemap = GeoData.GeoArray(nodemap, dims = lat_lon_dims)

    nodemap, paths
end

function path_to_array(path::Vector{Int},
                       nodemap::Matrix{Int})
    path_array = zeros(eltype(nodemap), size(nodemap))
    path_array[in.(nodemap, [path])] .= 1

    return path_array
end

function path_to_geoarray(path::Vector{Int},
                          nodemap::GeoData.GeoArray)
    path_array = path_to_array(path, nodemap.data[:, :, 1])
    lat_lon_dims = get_lat_lon_dims(weights)

    return GeoData.GeoArray(path_array, dims = lat_lon_dims)
end

# Convert a path to a vector of its coordinates. Provide a geotransform to
# get proper geographic coordinates
function path_to_points(path::Vector{Int},
                        nodemap::Matrix{Int};
                        geotransform::Vector{N} where N <: Real = [0.0, 1.0, 0.0, 0.0, 0.0, -1.0],
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

