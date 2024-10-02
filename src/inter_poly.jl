include("derivation.jl")

function eval_polynom(polynom::AbstractVector, point::Number)
    powers = point.^(0:length(polynom)-1)
    return polynom'powers
end

function root_polynom(points::T) where T <: Union{AbstractVector, AbstractSet}
    isempty(points) && return [1]
    poly = zeros(length(points) + 1)
    poly[end] = one(first(points))
    for (power, point) in enumerate(points)
        poly[end-power:end-1] -= point * poly[end-power+1:end] 
    end
    return poly
end

function lagrange_polynom(points::AbstractVector, vals::AbstractVector)
    length(points) == length(vals) || throw(ArgumentError("different x and y lengths"))
    isempty(points) && return []
    length(points) == 1 && return [vals[1]]
    points_set = Set(points)
    length(points_set) == length(points) || throw(ArgumentError("interpolation through same points"))
    inter_poly = zeros(length(points))
    for (point, val) in zip(points, vals)
        delete!(points_set, point)
        poly_without_point = root_polynom(points_set)
        inter_poly += val * poly_without_point / eval_polynom(poly_without_point, point) 
        push!(points_set, point)
    end
    return inter_poly
end

function newton_polynom(points::AbstractVector, vals::AbstractVector)
    length(points) == length(vals) || throw(ArgumentError("different x and y lengths"))
    isempty(points) && return []
    length(points) == 1 && return [vals[1]]
    points_set = Set(points)
    length(points_set) == length(points) || throw(ArgumentError("interpolation through same points"))
    inter_poly = zeros(length(points))
    inter_poly[begin] = vals[1]
    for power in 1:(length(points)-1)
        inter_poly[begin:begin+power] .+= root_polynom(points[begin:begin+power-1]) *
            divided_diff(points[begin:begin+power], vals[begin:begin+power])
    end
    return inter_poly
end