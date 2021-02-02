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

    cd_geoarray = GeoData.GeoArray(reshape(cd_array, size(nodemap)),
                                   dims = lat_lon_dims,
                                   missingval = -1)
    cd_geoarray
end

function sample_node_pairs(sample_weights::Matrix{T} where T <: Real, # Matrix of weights
                               nodemap::Matrix{Int}, # Matrix of node_ids to sample
                               n_pairs::Int)
    weight = Float64.(deepcopy(sample_weights))
    # Set weights (habitat) to 0 for non-nodes
    weight[nodemap .== 0] .= 0

    # Mask weights objects from habitat layer for use in sampling
    weight_weights = StatsBase.weights(weight)

    # Initialize object in which to store pairs of nodemap
    samples = Vector{Vector{Int}}()

    # Sample n_pairs * 2 points, sample apparently needs Float 64, so coerced above
    the_sample = sample(nodemap, weight_weights, n_pairs * 2, replace = true)

    for i in 1:(n_pairs)
        push!(samples, the_sample[(2 * i - 1):(2 * i)])
    end

    return samples
end

function sample_node_pairs(sample_weights::GeoData.GeoArray,
                               nodemap::GeoData.GeoArray,
                               n_pairs::Int)
    weight = deepcopy(sample_weights.data[:, :, 1])
    weight[weight .== sample_weights.missingval] .= 0
    weight[isnan.(weight)] .= 0

    samples = sample_node_pairs(weight,
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
                     connect_using_avg_cost::Bool = true,
                     parallel::Bool = true)
    @info "Constructing graphs"
    nodemap = construct_nodemap(cost_surface)
    graph = construct_weighted_graph(cost_surface,
                            nodemap,
                            no_data_val = no_data_val,
                            cost_layer_is_conductance = cost_layer_is_conductance,
                            connect_four_neighbors_only = connect_four_neighbors_only,
                            connect_using_avg_cost = connect_using_avg_cost)

    @info "Generating random path start and end points"
    node_pairs = sample_node_pairs(sample_weights,
                                       nodemap,
                                       n_paths)
    #### Identify the least cost paths
    # Initialize the paths
    @info "Computing least cost paths"
    paths = Vector{Vector{Int}}(undef, n_paths)
    p = Progress(n_paths;
                 dt = 0.25,
                 barlen = min(50, displaysize(stdout)[2] - length("Progress: 100%  Time: 00:00:00")))
    if parallel
        @threads for i in 1:n_paths
           lcp = least_cost_path(graph, node_pairs[i][1], node_pairs[i][2])
           paths[i] = lcp
           next!(p)
        end
    else
        for i in 1:n_paths
           lcp = least_cost_path(graph, node_pairs[i][1], node_pairs[i][2])
           paths[i] = lcp
           next!(p)
        end
    end

    return paths, nodemap
end

function random_lcps(cost_surface::GeoData.GeoArray,
                     sample_weights::GeoData.GeoArray,
                     n_paths::Int;
                     cost_layer_is_conductance::Bool = false,
                     connect_four_neighbors_only::Bool = false,
                     connect_using_avg_cost::Bool = true,
                     parallel::Bool = true)
    @info "Constructing graphs"
    nodemap = construct_nodemap(cost_surface)
    graph = construct_weighted_graph(cost_surface,
                            nodemap,
                            cost_layer_is_conductance = cost_layer_is_conductance,
                            connect_four_neighbors_only = connect_four_neighbors_only,
                            connect_using_avg_cost = connect_using_avg_cost)

    @info "Generating random path start and end points"
    node_pairs = sample_node_pairs(sample_weights,
                                       nodemap,
                                       n_paths)
    #### Identify the least cost paths
    # Initialize the paths
    @info "Computing least cost paths"
    paths = Vector{Vector{Int}}(undef, n_paths)
    p = Progress(n_paths;
                 dt = 0.25,
                 barlen = min(50, displaysize(stdout)[2] - length("Progress: 100%  Time: 00:00:00")))

    if parallel
        @threads for i in 1:n_paths
            lcp = least_cost_path(graph, node_pairs[i][1], node_pairs[i][2])
            paths[i] = lcp
            next!(p)
        end
    else
        for i in 1:n_paths
            lcp = least_cost_path(graph, node_pairs[i][1], node_pairs[i][2])
            paths[i] = lcp
            next!(p)
        end
    end

    return paths, nodemap
end

function path_to_array(path::Vector{Int},
                       nodemap::Matrix{Int})
    path_array = zeros(eltype(nodemap), size(nodemap))

    # nodemap needs to have been created by construct_nodemap or other
    # SpatialGraphs function for this to work properly
    coords = path_to_cartesian_coords(path, nodemap)
    path_array[coords] .= 1

    return path_array
end

function path_to_geoarray(path::Vector{Int},
                          nodemap::GeoData.GeoArray)
    path_array = reshape(path_to_array(path, nodemap.data[:, :, 1]),
                         size(nodemap.data))
    lat_lon_dims = get_lat_lon_dims(nodemap)

    return GeoData.GeoArray(path_array, dims = lat_lon_dims)
end

# Convert a path to a vector of its coordinates.
function path_to_cartesian_coords(path::Vector{Int},
                                  nodemap::Matrix{Int})
    cart_coords = get_cartesian_indices(nodemap)[path]

    return cart_coords
end

function path_to_linestring(path::Vector{Int},
                            nodemap::GeoData.GeoArray)
    # first get cartesian coordinates of the nodemap
    cart_coords = get_cartesian_indices(nodemap)[path]

    # Convert to geo coordinates
    lat_lon_dims = SpatialGraphs.get_lat_lon_dims(nodemap)
    row_dim = lat_lon_dims[1].val
    row_step = step(lat_lon_dims[1]) * 0.5
    col_dim = lat_lon_dims[2].val
    col_step = step(lat_lon_dims[2]) * 0.5

    # Robust to permuations of GeoArray
    lon_first = typeof(lat_lon_dims[1]) <: Lon
    if lon_first
        geo_coords = map(x -> (row_dim[x[1]] + row_step, col_dim[x[2]] - col_step), cart_coords)
    else
        geo_coords = map(x -> (col_dim[x[2]] - col_step, row_dim[x[1]] + row_step), cart_coords)
    end

    return createlinestring(geo_coords)
end

function paths_to_linestrings(paths::Vector{Vector{Int}},
                              nodemap::GeoData.GeoArray;
                              parallel = true)
    # Info for conversion to geo coordinates
    lat_lon_dims = SpatialGraphs.get_lat_lon_dims(nodemap)
    row_dim = lat_lon_dims[1].val
    row_step = step(lat_lon_dims[1]) * 0.5
    col_dim = lat_lon_dims[2].val
    col_step = step(lat_lon_dims[2]) * 0.5

    # get cartesian coordinates
    all_cart_coords = get_cartesian_indices(nodemap)
    geo_coords_all = Vector{Vector{Tuple{Float64, Float64}}}(undef, length(paths))

    lon_first = typeof(lat_lon_dims[1]) <: Lon

    if parallel
        @threads for i in 1:length(paths)
            cart_coords = all_cart_coords[paths[i]]
            if lon_first
                geo_coords = map(x -> (row_dim[x[1]] + row_step, col_dim[x[2]] - col_step), cart_coords)
                geo_coords_all[i] = geo_coords
            else
                geo_coords = map(x -> (col_dim[x[2]] - col_step, row_dim[x[1]] + row_step), cart_coords)
                geo_coords_all[i] = geo_coords
            end
        end
    else
        for i in 1:length(paths)
            cart_coords = all_cart_coords[paths[i]]
            if lon_first
                geo_coords = map(x -> (row_dim[x[1]] + row_step, col_dim[x[2]] - col_step), cart_coords)
                geo_coords_all[i] = geo_coords
            else
                geo_coords = map(x -> (col_dim[x[2]] - col_step, row_dim[x[1]] + row_step), cart_coords)
                geo_coords_all[i] = geo_coords
            end
        end
    end

    linestrings = createlinestring.(geo_coords_all)

    return linestrings
end
