include("../src/OVF.jl")
using .OVF

space_grid = range(
    start=-3,
    stop=3,
    length=10
)
potential(x) = x^2 / 2
potential_values = potential.(space_grid)
energy = minimum(potential_values)
Ψ_state = exp.(energy .- potential_values)
Ψ_state /= sqrt(Ψ_state'Ψ_state)

Δx = step(space_grid)
top = bot = fill(-1/2Δx^2, length(space_grid) - 1)

for _ in 1:10

end

println(energy)
using Plots
plot(space_grid, Ψ_state)
png("plots/task_11_wave_function.png")