module SpatialGraphs
using LightGraphs, SimpleWeightedGraphs, GeoData

include("graphs.jl")
export construct_nodemap, construct_graph

include("paths.jl")
export sample_lcp_node_pairs, compute_cd_graph, least_cost_path

end # module