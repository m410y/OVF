include("../src/OVF.jl")
using .OVF

grid = range(0, 2pi, 300_000)
j0_array = OVF.j0.(grid)
j1_array = OVF.j1.(grid)
j0_deriv = OVF.derivative(j0_array, grid)
task_check = j0_deriv .+ j1_array

using Plots
plot(grid, task_check)
png("plots/task_4_j0_diff_plus_j1.png")