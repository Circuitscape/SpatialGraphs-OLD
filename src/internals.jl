# averaging functions, to operate on resistances
# cond_* functions convert to conductance, calulate mean, then convert to resistance
# res_* functions directly operate on the resistances
cond_cardinal_avg(x, y) = 1 / ((1/x + 1/y) / 2)
cond_diagonal_avg(x, y) = 1 / ((1/x + 1/y) / (2 * √2))
res_cardinal_avg(x, y) = (x + y) / 2
res_diagonal_avg(x, y) = ((x + y) * √2) / 2


function get_lat_lon_dims(A::GeoData.GeoArray; layer_name::String = "input")
    dim_names = String.(name(dims(A)))
    lat_first = findall(dim_names .== "Latitude")[1] < findall(dim_names .== "Longitude")[1]
    first_dim_type = lat_first ? Lat : Lon
    second_dim_type = lat_first ? Lon : Lat

    (dims(resistance, first_dim_type), dims(resistance, second_dim_type))
end