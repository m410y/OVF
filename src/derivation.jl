function derivative(arr::AbstractVector, grid::AbstractRange)
    length(arr) == length(grid) || throw(ArgumentError("different length"))
    left = -1.5arr[begin] + 2arr[begin+1] - 0.5arr[begin+2]
    right = 0.5arr[end-2] - 2arr[end-1] + 1.5arr[end]
    middle = 0.5(arr[begin+2:end] - arr[begin:end-2])
    return vcat(left, middle, right) ./ step(grid)
end

function divided_diff(points::AbstractVector, vals::AbstractVector)
    length(points) == length(vals) || throw(ArgumentError("different length"))
    for depth in 1:(length(points) - 1)
        vals = (vals[begin+1:end] - vals[begin:end-1]) ./
            (points[begin+depth:end] - points[begin:end-depth])
    end
    return first(vals)
end