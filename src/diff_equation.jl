include("root_find.jl")
include("linear_algebra.jl")

function euler_method(f::Function, grid::AbstractRange, start_value::T) where {T<:Number}
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

function euler_method(
    f::Function,
    grid::AbstractRange,
    start_value::AbstractVector{T},
) where {T<:Number}
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

function rk2_method(
    f::Function,
    grid::AbstractRange,
    start_value::T;
    α = 1,
) where {T<:Number}
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

function rk2_method(
    f::Function,
    grid::AbstractRange,
    start_value::AbstractVector{T};
    α = 1,
) where {T<:Number}
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

function rk4_method(f::Function, grid::AbstractRange, start_value::T) where {T<:Number}
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

function rk4_method(
    f::Function,
    grid::AbstractRange,
    start_value::AbstractVector{T},
) where {T<:Number}
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

function euler_implicit_method(
    f::Function,
    J::Function,
    grid::AbstractRange,
    start_value::T;
    α = 1,
    iterations = 1,
) where {T<:Number}
    solution = Vector{T}(undef, length(grid))
    h = step(grid)
    y = start_value
    solution[begin] = y
    for (step_n, x) in enumerate(grid[begin:end-1])
        y += newton_root(
            t -> t - (1 - α) * h * f(x, y) - α * h * f(x + h, y + t),
            t -> 1 - α * h * J(x + h, y + t),
            zero(T),
            iter_max = iterations,
        )
        solution[step_n+1] = y
    end
    return solution
end

function euler_implicit_method(
    f::Function,
    J::Function,
    grid::AbstractRange,
    start_value::AbstractVector{T};
    α = 1,
    iterations = 1,
) where {T<:Number}
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
            iter_max = iterations,
        )
        solution[step_n+1, :] = y
    end
    return solution
end

function poisson_1D(
    f::Function,
    grid::AbstractRange,
    left_border::Tuple,
    right_border::Tuple,
)
    N = length(grid)
    A_l, B_l, C_l = left_border
    A_r, B_r, C_r = right_border
    h = step(grid)
    bot = [fill(1 / h^2, N - 2); -B_r / h]
    diag = [A_l - B_l / h; fill(-2 / h^2, N - 2); A_r + B_r / h]
    top = [B_l / h; fill(1 / h^2, N - 2)]
    vals = [C_l; f.(grid[begin+1:end-1]); C_r]
    return thomas_algorithm!(bot, diag, top, vals)
end

function conductivity_1D(
    f::Function,
    space_grid::AbstractRange,
    time_grid::AbstractRange,
    left_border::Tuple,
    right_border::Tuple,
    init_vals::AbstractVector{T};
    α = 1,
) where {T<:Real}
    N_t, N_x = length.((time_grid, space_grid))
    length(init_vals) == N_x || throw(ArgumentError("wrong initial data size"))
    A_l, B_l, C_l = left_border
    A_r, B_r, C_r = right_border
    Δx = step(space_grid)
    Δt = step(time_grid)
    u = Matrix{T}(undef, N_t, N_x)
    u[begin, :] = init_vals
    bot = [fill(-α / 2Δx^2, N_x - 2); -B_r / Δx]
    diag = [
        A_l - B_l / Δx
        fill(1 / Δt + α / Δx^2, N_x - 2)
        A_r + B_r / Δx
    ]
    top = [B_l / Δx; fill(-α / 2Δx^2, N_x - 2)]
    vals = [C_l; Vector{T}(undef, N_x - 2); C_r]
    for (t_n, t) in enumerate(time_grid[begin:end-1])
        u_now = u[t_n, :]
        vals[begin+1:end-1] =
            map(x -> f(x, t + Δt / 2), space_grid[begin+1:end-1]) +
            u_now[begin:end-2] * α / 2Δx^2 +
            u_now[begin+1:end-1] * (1 / Δt - α / Δx^2) +
            u_now[begin+2:end] * α / 2Δx^2
        u[t_n+1, :] = thomas_algorithm(bot, diag, top, vals)
    end
    return u
end

function schrödinger_stationary_1D(
    U::Function,
    grid::AbstractRange;
    E_0 = minimum(U.(grid)),
    Ψ_0 = exp.(E_0 .- U.(grid)),
    error = eps(E_0),
    iter_max = 9999,
)
    N = length(grid)
    length(Ψ_0) == N || throw(ArgumentError("wrong Ψ and grid length"))
    h = step(grid)
    U_vals = U.(grid)
    Ψ = Ψ_0 / sqrt(trapezoid_integral(Ψ_0 .^ 2, h))
    E_old, E = -Inf, E_0
    top = bot = fill(-1 / 2h^2, N - 1)
    diag = fill(1 / h^2, N) .+ U_vals
    iter = 0
    while abs(E_old - E) > error && iter < iter_max
        Ψ = thomas_algorithm(bot, diag, top, E * Ψ)
        E_old, E = E, sqrt(trapezoid_integral(Ψ .^ 2, h))
        Ψ /= E
        iter += 1
    end
    return E, Ψ
end
