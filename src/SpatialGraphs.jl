module SpatialGraphs
using LightGraphs, SimpleWeightedGraphs, StatsBase,
      Base.Threads, ProgressMeter, GeoData
using ArchGDAL: createlinestring

include("internals.jl")

include("graphs.jl")
export construct_nodemap, construct_weighted_graph

include("paths.jl")
export sample_node_pairs, cost_distance, least_cost_path,
       path_to_cartesian_coords, random_lcps, path_to_linestring,
       paths_to_linestrings, path_to_array, path_to_geoarray

end
