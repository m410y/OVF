function derivative(vals::AbstractVector{T}, h::T) where {T<:Number}
    left = -3vals[begin] + 4vals[begin+1] - vals[begin+2]
    right = vals[end-2] - 4vals[end-1] + 3vals[end]
    middle = vals[begin+2:end] - vals[begin:end-2]
    return [left; middle; right] ./ 2h
end

function divided_diff(vals::AbstractVector{T}, points::AbstractVector{T}) where {T<:Number}
    N = length(vals)
    length(points) == N || throw(ArgumentError("different length"))
    for depth = 1:(N-1)
        vals =
            (vals[begin+1:end] - vals[begin:end-1]) ./
            (points[begin+depth:end] - points[begin:end-depth])
    end
    return first(vals)
end
