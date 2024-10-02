include("root_find.jl")
include("linear_algebra.jl")

function euler_method(f::Function, grid::AbstractRange, start_value::T) where T <: Number
    solution = Vector{T}(undef, length(grid))
    h = step(grid)
    y = start_value
    solution[begin] = y
    for (step_n, x) in enumerate(grid[begin:end-1])
        y += f(x, y) * h
        solution[step_n+1] = y
    end
    return solution
end

function euler_method(f::Function, grid::AbstractRange, start_value::AbstractVector{T}) where T <: Number
    solution = Matrix{T}(undef, length(grid), length(start_value))
    h = step(grid)
    y = start_value
    solution[begin, :] = y
    for (step_n, x) in enumerate(grid[begin:end-1])
        y += f(x, y) * h
        solution[step_n+1, :] = y
    end
    return solution
end

function rk2_method(f::Function, grid::AbstractRange, start_value::T; α=1) where T <: Number
    solution = Vector{T}(undef, length(grid))
    h = step(grid)
    y = start_value
    solution[begin] = y
    for (step_n, x) in enumerate(grid[begin:end-1])
        k_1 = f(x, y)
        k_2 = f(x + h / 2α, y + h / 2α * k_1)
        y += ((1 - α) * k_1 + α * k_2) * h
        solution[step_n+1] = y
    end
    return solution
end

function rk2_method(f::Function, grid::AbstractRange, start_value::AbstractVector{T}; α=1) where T <: Number
    solution = Matrix{T}(undef, length(grid), length(start_value))
    h = step(grid)
    y = start_value
    solution[begin, :] = y
    for (step_n, x) in enumerate(grid[begin:end-1])
        k_1 = f(x, y)
        k_2 = f(x + h / 2α, y + h / 2α * k_1)
        y += ((1 - α) * k_1 + α * k_2) * h
        solution[step_n+1, :] = y
    end
    return solution
end

function rk4_method(f::Function, grid::AbstractRange, start_value::T) where T <: Number
    solution = Vector{T}(undef, length(grid))
    h = step(grid)
    y = start_value
    solution[begin] = y
    for (step_n, x) in enumerate(grid[begin:end-1])
        k_1 = f(x, y)
        k_2 = f(x + h / 2, y + h / 2 * k_1)
        k_3 = f(x + h / 2, y + h / 2 * k_2)
        k_4 = f(x + h, y + h * k_3)
        y += (k_1 + 2k_2 + 2k_3 + k_4) * h / 6
        solution[step_n+1] = y
    end
    return solution
end

function rk4_method(f::Function, grid::AbstractRange, start_value::AbstractVector{T}) where T <: Number
    solution = Matrix{T}(undef, length(grid), length(start_value))
    h = step(grid)
    y = start_value
    solution[begin, :] = y
    for (step_n, x) in enumerate(grid[begin:end-1])
        k_1 = f(x, y)
        k_2 = f(x + h / 2, y + h / 2 * k_1)
        k_3 = f(x + h / 2, y + h / 2 * k_2)
        k_4 = f(x + h, y + h * k_3)
        y += (k_1 + 2k_2 + 2k_3 + k_4) * h / 6
        solution[step_n+1, :] = y
    end
    return solution
end

function euler_implicit_method(f::Function, J::Function, grid::AbstractRange, start_value::T; α=1, iterations=1) where T <: Number
    solution = Vector{T}(undef, length(grid))
    h = step(grid)
    y = start_value
    solution[begin] = y
    for (step_n, x) in enumerate(grid[begin:end-1])
        y += newton_root(
            t -> t - (1 - α) * h * f(x, y) - α * h * f(x + h, y + t),
            t -> 1 - α * h * J(x + h, y + t),
            zero(y),
            max_iter=iterations
        )
        solution[step_n+1] = y
    end
    return solution
end

function euler_implicit_method(f::Function, J::Function, grid::AbstractRange, start_value::AbstractVector{T}; α=1, iterations=1) where T <: Number
    solution = Matrix{T}(undef, length(grid), length(start_value))
    h = step(grid)
    y = start_value
    solution[begin, :] = y
    I = identity_matrix(length(y))
    for (step_n, x) in enumerate(grid[begin:end-1])
        y += newton_root(
            t -> t - (1 - α) * h * f(x, y) - α * h * f(x + h, y + t),
            t -> I - α * h * J(x + h, y + t),
            zero(y),
            max_iter=iterations
        )
        solution[step_n+1, :] = y
    end
    return solution
end

function poisson_1D_bounded(f::Function, grid::AbstractRange, left_border::Tuple, right_border::Tuple)
    A_left, B_left, C_left = left_border
    A_right, B_right, C_right = right_border
    h = step(grid)
    one_line =  fill(1, length(grid) - 2)
    solution = thomas_algorithm!(
        [one_line/h^2; -B_right/h],
        [A_left - B_left/h; -2one_line/h^2; A_right + B_right/h],
        [B_left/h; one_line/h^2],
        [C_left; f.(grid[begin+1:end-1]); C_right]
    )
    return solution
end