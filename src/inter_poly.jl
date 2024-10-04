include("derivation.jl")

function eval_polynom(polynom::AbstractVector{T}, point::Number) where {T<:Number}
    powers = point .^ (0:length(polynom)-1)
    return polynom'powers
end

function root_polynom(points::Union{AbstractVector{T},AbstractSet{T}}) where {T<:Number}
    isempty(points) && return [1]
    poly = zeros(length(points) + 1)
    poly[end] = one(first(points))
    for (power, point) in enumerate(points)
        poly[end-power:end-1] -= point * poly[end-power+1:end]
    end
    return poly
end

function lagrange_polynom(
    vals::AbstractVector{T},
    points::AbstractVector{T},
) where {T<:Number}
    N = length(vals)
    length(points) == N || throw(ArgumentError("different x and y lengths"))
    N == 0 && return []
    N == 1 && return [first(vals)]
    points_set = Set(points)
    length(points_set) == N || throw(ArgumentError("interpolation through same points"))
    inter_poly = zeros(N)
    for (point, val) in zip(points, vals)
        delete!(points_set, point)
        poly_without_point = root_polynom(points_set)
        inter_poly += val * poly_without_point / eval_polynom(poly_without_point, point)
        push!(points_set, point)
    end
    return inter_poly
end

function newton_polynom(
    vals::AbstractVector{T},
    points::AbstractVector{T},
) where {T<:Number}
    N = length(vals)
    length(points) == N || throw(ArgumentError("different x and y lengths"))
    N == 0 && return []
    N == 1 && return [first(vals)]
    points_set = Set(points)
    length(points_set) == N || throw(ArgumentError("interpolation through same points"))
    inter_poly = zeros(N)
    inter_poly[begin] = vals[begin]
    for power = 1:(N-1)
        inter_poly[begin:begin+power] .+=
            root_polynom(points[begin:begin+power-1]) *
            divided_diff(vals[begin:begin+power], points[begin:begin+power])
    end
    return inter_poly
end
