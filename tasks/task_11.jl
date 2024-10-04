include("../src/OVF.jl")
using .OVF

space_grid = range(start = -3, stop = 3, length = 10)

energy, wave_function = OVF.schrödinger_stationary_1D(x -> x^2 / 2, space_grid)

println(energy)
using Plots
plot(space_grid, wave_function)
png("plots/task_11_wave_function.png")
