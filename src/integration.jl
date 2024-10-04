function middle_integral(vals::AbstractVector{T}, h::T) where {T<:Number}
    return sum(vals) * h
end

function trapezoid_integral(vals::AbstractVector{T}, h::T) where {T<:Number}
    return (2sum(vals) - vals[begin] - vals[end]) * h / 2
end

function simpson_integral(vals::AbstractVector{T}, h::T) where {T<:Number}
    isodd(length(vals)) || throw(ArgumentError("even grid"))
    return (
        2sum(vals[begin:2:end]) + 4sum(vals[begin+1:2:end-1]) - vals[begin] - vals[end]
    ) * h / 3
end

function left_integral(f::Function, grid::AbstractRange)
    return middle_integral(f.(grid[begin:end-1]), step(grid))
end

function right_integral(f::Function, grid::AbstractRange)
    return middle_integral(f.(grid[begin+1:end]), step(grid))
end

function middle_integral(f::Function, grid::AbstractRange)
    middle_grid = range(
        start = grid[begin] + step(grid) / 2,
        step = step(grid),
        length = length(grid) - 1,
    )
    return middle_integral(f.(middle_grid), step(middle_grid))
end

function trapezoid_integral(f::Function, grid::AbstractRange)
    return trapezoid_integral(f.(grid), step(grid))
end

function simpson_integral(f::Function, grid::AbstractRange)
    odd_grid = range(
        start = grid[begin],
        stop = grid[end],
        length = length(grid) + (length(grid) + 1) % 2,
    )
    return simpson_integral(f.(odd_grid), step(odd_grid))
end
