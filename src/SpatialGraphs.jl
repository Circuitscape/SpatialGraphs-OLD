module SpatialGraphs
using LightGraphs, SimpleWeightedGraphs, GeoData, Base.Threads

include("internals.jl")

include("graphs.jl")
export construct_nodemap, construct_graph

include("paths.jl")
export sample_lcp_node_pairs, cost_distance, least_cost_path, path_to_array,
       path_to_points, random_lcps, cd_to_array, cd_to_geoarray

end
