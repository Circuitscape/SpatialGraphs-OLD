# TODO add a function that takes output from this and converts it to a GeoArray or Matrix
function cost_distance(g::AbstractSimpleWeightedGraph,
                       nodemap::Union{Matrix{Int}, GeoData.GeoArray},
                       node_ids::Union{Int, Array{Int, 1}}; # will report lowest CD to get to any of these, Number or Vector
                       cost_threshold::Real = Inf,
                       dist_fun::Function = dijkstra_shortest_paths)
    pathstate = dist_fun(g, node_ids)

    # Return type based on nodemap type
    if isa(nodemap, Matrix{Int})
        return pathstate_to_array(pathstate, nodemap,
                                  cost_threshold = cost_threshold)
    else
        return pathstate_to_geoarray(pathstate, nodemap,
                                     cost_threshold = cost_threshold)
    end
end

function pathstate_to_array(cdist::LightGraphs.AbstractPathState,
                            nodemap::Matrix{Int};
                            cost_threshold::Real = Inf)
    cd_array = fill(-1., size(nodemap))
    cd_array[nodemap .!= 0] .= cdist.dists[nodemap[nodemap.!= 0]]

    cd_array[cd_array .> cost_threshold] .= -1

    cd_array
end

function pathstate_to_geoarray(cdist::LightGraphs.AbstractPathState,
                               nodemap::GeoData.GeoArray;
                               cost_threshold::Real = Inf)
    cd_array = pathstate_to_array(cdist,
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
    weight = float.(deepcopy(sample_weights))
    # Set weights (habitat) to 0 for non-nodes
    weight[nodemap .== 0] .= 0

    # Mask weights objects from habitat layer for use in sampling
    weight_weights = StatsBase.weights(weight)

    # Initialize object in which to store pairs of nodemap
    samples = Vector{Vector{Int}}()

    # Sample n_pairs * 2 points
    the_sample = sample(nodemap, weight_weights, n_pairs * 2, replace = true)

    for i in 1:(n_pairs)
        push!(samples, the_sample[(2 * i - 1):(2 * i)])
    end

    return samples
end

function sample_lcp_node_pairs(sample_weights::GeoData.GeoArray,
                               nodemap::GeoData.GeoArray,
                               n_pairs::Int)
    weight = deepcopy(sample_weights.data[:, :, 1])
    weight[weight .== sample_weights.missingval] .= 0
    weight[isnan.(weight)] .= 0
    
    samples = sample_lcp_node_pairs(weight,
                                    nodemap.data[:, :, 1],
                                    n_pairs)

    return samples
end

function least_cost_path(g::AbstractGraph,
                         start::Union{Int, Vector{T} where T <: Int},
                         destination::Int;
                         dist_fun::Function = dijkstra_shortest_paths)
    pathstate = dist_fun(g, start)
    return enumerate_paths(pathstate, destination)
end

# Returns a nodemap to map node_ids to points in space
# and a vector of vectors (each vector is a path)
function random_lcps(cost_surface::Matrix{T} where T <: Real,
                     sample_weights::Matrix{T} where T <: Real,
                     n_paths::Int;
                     no_data_val::Union{Nothing, Real} = nothing,
                     cost_layer_is_conductance::Bool = false,
                     connect_four_neighbors_only::Bool = false,
                     connect_using_avg_resistance::Bool = false,
                     parallel::Bool = true)
    nodemap = construct_nodemap(cost_surface)
    graph = construct_graph(cost_surface,
                            nodemap,
                            no_data_val = no_data_val,
                            cost_layer_is_conductance = cost_layer_is_conductance,
                            connect_four_neighbors_only = connect_four_neighbors_only,
                            connect_using_avg_resistance = connect_using_avg_resistance)

    node_pairs = sample_lcp_node_pairs(sample_weights,
                                       nodemap,
                                       n_paths)
    #### Identify the least cost paths
    # Initialize the paths
    paths = Vector{Vector{Int}}(undef, n_paths)

    if parallel
        @threads for i in 1:n_paths
           lcp = least_cost_path(graph, node_pairs[i][1], node_pairs[i][2])
           paths[i] = lcp
        end
    else
        for i in 1:n_paths
           lcp = least_cost_path(graph, node_pairs[i][1], node_pairs[i][2])
           paths[i] = lcp
        end
    end

    return nodemap, paths
end


function random_lcps(cost_surface::GeoData.GeoArray,
                     sample_weights::GeoData.GeoArray,
                     n_paths::Int;
                     cost_layer_is_conductance::Bool = false,
                     connect_four_neighbors_only::Bool = false,
                     connect_using_avg_resistance::Bool = false,
                     parallel::Bool = true)
    nodemap, paths = random_lcps(cost_surface.data[:, :, 1],
                                 sample_weights.data[:, :, 1],
                                 n_paths;
                                 no_data_val = cost_surface.missingval,
                                 cost_layer_is_conductance = cost_layer_is_conductance,
                                 connect_four_neighbors_only = connect_four_neighbors_only,
                                 connect_using_avg_resistance = connect_using_avg_resistance,
                                 parallel = parallel)
    # Convert nodemap to GeoArray
    lat_lon_dims = get_lat_lon_dims(cost_surface)
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
    lat_lon_dims = get_lat_lon_dims(nodemap)

    return GeoData.GeoArray(path_array, dims = lat_lon_dims)
end

# Convert a path to a vector of its coordinates. Provide a geotransform to
# get proper geographic coordinates
function path_to_cartesian_coords(path::Vector{Int},
                                  nodemap::Matrix{Int};
                                  parallel::Bool = true)
    cart_coords = Vector{CartesianIndex{2}}(undef, length(path))

    if parallel
        @threads for i in 1:length(path)
            cart_coord = findall(nodemap .== path[i])
            cart_coords[i] = cart_coord[1]
        end
    else
        for i in 1:length(path)
            cart_coord = findall(nodemap .== path[i])
            cart_coords[i] = cart_coord[1]
        end
    end

    return cart_coords
end

function path_to_linestring(path::Vector{Int},
                            nodemap::GeoData.GeoArray;
                            parallel = true)
    # first get cartesian coordinates of the nodemap
    cart_coords = path_to_cartesian_coords(path,
                                           nodemap.data[:, :, 1],
                                           parallel = true)
    # Convert to geo coordinates
    # Convert to geo coordinates
    lat_lon_dims = SpatialGraphs.get_lat_lon_dims(nodemap)
    row_dim = lat_lon_dims[1].val
    row_step = step(lat_lon_dims[1]) * 0.5
    col_dim = lat_lon_dims[2].val
    col_step = step(lat_lon_dims[2]) * 0.5

    # Robust to permuations of GeoArray
    if typeof(lat_lon_dims[1]) <: Lon
        geo_coords = map(x -> (row_dim[x[1]] + row_step, col_dim[x[2]] - col_step), cart_coords)
    else
        geo_coords = map(x -> (col_dim[x[2]] - col_step, row_dim[x[1]] + row_step), cart_coords)
    end

    return createlinestring(geo_coords)
end
