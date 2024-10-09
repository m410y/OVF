include("../src/OVF.jl")
using .OVF
using Plots, Measures

L = 1
space_grid = range(start = 0, stop = L, length = 20)
time_grid = range(start = 0, stop = 1, length = 100)
left_border = (1, 0, 0)
right_border = (1, 0, 0)
initial_values = sin.(pi / L * space_grid)

solution = OVF.conductivity_1D(
    (x, t) -> 0.0,
    space_grid,
    time_grid,
    left_border,
    right_border,
    initial_values,
    α = 1,
)

solution_true = exp.(-time_grid .* (pi / L)^2) * initial_values'
heatmap(space_grid, time_grid, solution - solution_true, right_margin = 10.0mm)
png("plots/task_10_solution_difference_heatmap.png")
