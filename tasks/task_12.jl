include("../src/OVF.jl")
using .OVF
using Plots

time_grid = range(start = 0, stop = 2pi, length = 64)
frequency_grid = range(
    start = 0,
    step = 2pi / (time_grid[end] - time_grid[begin]),
    length = Int(ceil(length(time_grid) / 2)),
)
f(t) = sin(5.1t) + 0.002sin(25.5t)
values = f.(time_grid)
fourier = OVF.DFT(values)
fourier_cropped = fourier[begin:begin+length(frequency_grid)-1]

plot(time_grid, values)
png("plots/task_12_initial.png")

plot(frequency_grid, abs2.(fourier_cropped))
png("plots/task_12_fourier.png")
