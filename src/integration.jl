function left_integral(f::Function, grid::AbstractRange)
    return sum(f, grid[begin:end-1]) * step(grid)
end

function right_integral(f::Function, grid::AbstractRange)
    return sum(f, grid[begin+1:end]) * step(grid)
end

function trapezoid_integral(f::Function, grid::AbstractRange)
    left = f(grid[begin])
    right = f(grid[end])
    int_sum = sum(f, grid[begin+1:end-1])
    return (left + 2int_sum + right) * step(grid) / 2
end

function middle_integral(f::Function, grid::AbstractRange)
    middle_grid = range(
        start=grid[begin] + step(grid)/2,
        step=step(grid),
        length=length(grid) - 1
    )
    return sum(f, middle_grid) * step(grid)
end

function simpson_integral(f::Function, grid::AbstractRange)
    odd_grid = range(
        start=grid[begin],
        stop=grid[end], 
        length=length(grid) + (length(grid) + 1) % 2
    )
    left = f(odd_grid[begin])
    right = f(odd_grid[end])
    int_sum_1 = sum(f, odd_grid[begin+1:2:end-1])
    int_sum_2 = sum(f, odd_grid[begin:2:end])
    return (-left + 4int_sum_1 + 2int_sum_2 - right) * step(odd_grid) / 3
end