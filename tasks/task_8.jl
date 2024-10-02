include("../src/OVF.jl")
using .OVF

diff_eq_grid = 0:2e-2:2
A = [
    998 1998;
    -999 -1999
]
start_value = [1.0, 0.0]
solution_true = Matrix{Float64}(undef, length(diff_eq_grid), 2)
solution_true[begin, :] = start_value
for (step_n, x) in enumerate(diff_eq_grid[begin:end-1])
    solution_true[begin+step_n, :] = exp(A * x) * start_value
end
solution_explicit = OVF.euler_method((x, y) -> A * y, diff_eq_grid, start_value)
solution_half = OVF.euler_implicit_method((x, y) -> A * y, (x, y) -> A, diff_eq_grid, start_value, α=1/2)
solution_implicit = OVF.euler_implicit_method((x, y) -> A * y, (x, y) -> A, diff_eq_grid, start_value, α=1)

using Plots
plot(xlims=(-0.1, 2.1), ylims=(-1.1, 0.1))
plot!(solution_true[:, 1], solution_true[:, 2])
plot!(solution_explicit[:, 1], solution_explicit[:, 2], marker=:circ)
plot!(solution_half[:, 1], solution_half[:, 2], marker=:circ)
plot!(solution_implicit[:, 1], solution_implicit[:, 2], marker=:circ)
png("plots/task_8_phase_solution.png")