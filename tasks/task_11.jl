include("../src/OVF.jl")
using .OVF
using Plots

space_grid = range(start = -5, stop = 5, length = 1000)
U(x) = x^2 / 2
let E = 0.5
    Ψ(x) = exp(-x^2 / 2) / pi^0.25
    global energy_true = E
    global wave_function_true = Ψ
end

energy, wave_function = OVF.schrödinger_stationary_1D(U, space_grid)

println(energy - energy_true)
plot(space_grid, wave_function - wave_function_true.(space_grid))
png("plots/task_11_wave_function.png")
