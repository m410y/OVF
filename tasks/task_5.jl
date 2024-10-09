include("../src/OVF.jl")
using .OVF

inter_grid = range(start = 0, stop = 10, length = 21)
values = OVF.j0.(inter_grid)
inter_poly = OVF.newton_polynom(values, inter_grid)

using Plots
plot(x -> OVF.eval_polynom(inter_poly, x) - OVF.j0(x), xlims = (0, 10), legend = false)
png("plots/task_5_j0_interpolation_difference.png")
