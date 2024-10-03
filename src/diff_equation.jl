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

function poisson_1D(f::Function, grid::AbstractRange, left_border::Tuple, right_border::Tuple)
    A_left, B_left, C_left = left_border
    A_right, B_right, C_right = right_border
    h = step(grid)
    bottom_line = [fill(h^-2, length(grid) - 2); -B_right/h]
    diagonal_line = [A_left - B_left/h; fill(-2h^-2, length(grid) - 2); A_right + B_right/h]
    top_line = [B_left/h; fill(h^-2, length(grid) - 2)]
    values = [C_left; f.(grid[begin+1:end-1]); C_right]
    solution = thomas_algorithm!(bottom_line, diagonal_line, top_line, values)
    return solution
end

function conductivity_1D(f::Function, space_grid::AbstractRange, time_grid::AbstractRange, left_border::Tuple, right_border::Tuple, initial_values::AbstractVector{T}; α=1) where T <: Real
    length(initial_values) == length(space_grid) || throw(ArgumentError("wrong initial data size"))
    A_left, B_left, C_left = left_border
    A_right, B_right, C_right = right_border
    Δx = step(space_grid)
    Δt = step(time_grid)
    isapprox(A_left * initial_values[begin] + B_left * (initial_values[begin+1] - initial_values[begin]) / Δx, C_left, atol=eps(T)) || throw(ArgumentError("left border condition not satisfied"))
    isapprox(A_right * initial_values[end] + B_right * (initial_values[end] - initial_values[end-1]) / Δx, C_right, atol=eps(T)) || throw(ArgumentError("right border condition not satisfied"))
    solution = Matrix{T}(undef, length(time_grid), length(space_grid))
    solution[begin, :] = initial_values
    bottom_line = [fill(-α/2Δx^2, length(space_grid) - 2); -B_right/Δx]
    diagonal_line = [A_left - B_left/Δx; fill(1/Δt + α/Δx^2, length(space_grid) - 2); A_right + B_right/Δx]
    top_line = [B_left/Δx; fill(-α/2Δx^2, length(space_grid) - 2)]
    values = [C_left; Vector{T}(undef, length(space_grid) - 2); C_right]
    for (t_n, t) in enumerate(time_grid[begin:end-1])
        values[begin+1:end-1] = map(x -> f(x, t + Δt/2), space_grid[begin+1:end-1])
        values[begin+1:end-1] += solution[t_n, begin:end-2] * α/2Δx^2
        values[begin+1:end-1] += solution[t_n, begin+1:end-1] * (1/Δt - α/Δx^2)
        values[begin+1:end-1] += solutionu[t_n, begin+2:end] * α/2Δx^2
        solution[t_n+1, :] = thomas_algorithm(bottom_line, diagonal_line, top_line, values)
    end
    return solution
end

function schrödinger_stationary_1D(U::Function, grid::AbstractRange; initial_energy=minimum(U.(grid)), initial_Ψ=exp.(initial_energy .- U.(grid)), error=eps(initial_energy), max_iter=9999)
    length(initial_Ψ) == length(grid) || throw(ArgumentError("wrong Ψ and grid length"))
    Ψ = initial_Ψ / sqrt(initial_Ψ'initial_Ψ)
    old_E = initial_energy
    top_line = bottom_line = fill(-1/2Δx^2, length(grid) - 2)
    diagonal_line = fill(1/Δx^2, length(grid) - 1)
end