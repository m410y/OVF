include("linear_algebra.jl")
include("integration.jl")

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
)
    N = length(grid)
    length(Ψ_0) == N || throw(ArgumentError("wrong Ψ and grid length"))
    h = step(grid)
    top = bot = fill(-1 / 2h^2, N - 1)
    diag = fill(1 / h^2, N) + U.(grid) .- E_0
    E, Ψ = inverse_tridiag_method(bot, diag, top, start = Ψ_0)
    Ψ /= sign(Ψ[end]) * sqrt(trapezoid_integral(Ψ .^ 2, h))
    return E + E_0, Ψ
end
