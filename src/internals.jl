# averaging functions, to operate on resistances
# cond_* functions convert to conductance, calulate mean, then convert to resistance
# res_* functions directly operate on the resistances
cond_cardinal_avg(x, y) = 1 / ((1/x + 1/y) / 2)
cond_diagonal_avg(x, y) = 1 / ((1/x + 1/y) / (2 * √2))
res_cardinal_avg(x, y) = (x + y) / 2
res_diagonal_avg(x, y) = ((x + y) * √2) / 2


function get_lat_lon_dims(A::GeoData.GeoArray)
    dim_names = String.(name(dims(A)))
    lat_first = findall(dim_names .== "Latitude")[1] < findall(dim_names .== "Longitude")[1]
    first_dim_type = lat_first ? Lat : Lon
    second_dim_type = lat_first ? Lon : Lat

    (dims(A, first_dim_type), dims(A, second_dim_type), Band(1:1, Categorical(order = Ordered())))
end

function nonunique!(x::Vector{Int})
    sort!(x)
    duplicatedvector = Int[]
    for i=2:length(x)
        if (isequal(x[i], x[i-1]) && (length(duplicatedvector) == 0 || !isequal(duplicatedvector[end], x[i])))
            push!(duplicatedvector,x[i])
        end
    end
    unique(duplicatedvector)
end

function edge_min_costs!(sources::Vector{Int},
                         destinations::Vector{Int},
                         node_weights::Vector{T} where T <: Real;
                         duplicated_node_ids::Vector{Int})
    ## Drop any self connections
    non_self_connected = sources .!= destinations
    sources = sources[non_self_connected]
    destinations = destinations[non_self_connected]
    node_weights = node_weights[non_self_connected]

    ## Sort edges (element-wise) so duplicates can be found
    sources = min.(sources, destinations)
    destinations = max.(sources, destinations)

    # Identify and separate out edges that can not have duplicates vs those that
    # could
    bad_edge_indices = in.(sources, [duplicated_node_ids]) .| in.(destinations, [duplicated_node_ids])
    good_sources = sources[(!).(bad_edge_indices)]
    good_dests = destinations[(!).(bad_edge_indices)]
    good_edge_weights = node_weights[(!).(bad_edge_indices)]

    dup_sources = sources[bad_edge_indices]
    dup_dests = destinations[bad_edge_indices]
    dup_edge_weights = node_weights[bad_edge_indices]

    # Initialize objects to store info
    idx = Dict{Tuple{Int,Int}, Int}((dup_sources[1], dup_dests[1]) => 1)
    weight_mins = [dup_edge_weights[1]]
    unique_sources = [dup_sources[1]]
    unique_dests = [dup_dests[1]]

    # For loop to identify the minimum cost posible to get from source to dest,
    # for each unique edge source-destination pair.
    Ct = 2
    @inbounds for i in 2:length(unique_sources)
        current_edge = (dup_sources[i], dup_dests[i])
        if haskey(idx, current_edge) # Has this edge already been found?
            # if yes, check if the current weight is lower than the one already
            # found, and update it if it is.
            weight_mins[idx[current_edge]] = min(weight_mins[idx[current_edge]], dup_edge_weights[i])
        else # Otherwise, push! the newly observed source, destination, and weight to the vectors storing them
            idx[current_edge] = Ct # Keep track of the index in weight_mins,
                                   # unique_sources, and unique_dests
                                   # corresponding to current_edge. It is
                                   # referenced later when weight_mins needs
                                   # a value overwritten because a new minimum
                                   # was found
            push!(weight_mins, dup_edge_weights[i])
            push!(unique_sources, dup_sources[i])
            push!(unique_dests, dup_dests[i])
            Ct += 1
        end
    end

    sources = vcat(good_sources, unique_sources)
    destinations = vcat(good_dests, unique_dests)
    node_weights = vcat(good_edge_weights, weight_mins)
end
