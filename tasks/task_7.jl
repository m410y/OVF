include("../src/OVF.jl")
using .OVF

diff_eq_grid = 0:1e-2:1
a, b, c, d = 10, 2, 2, 10
diff_eq_function = y -> [
    a * y[1] - b * y[1] * y[2],
    c * y[1] * y[2] - d * y[2]
]
start_value = [1.0, 1.0]
solution = OVF.rk2_method((x, y) -> diff_eq_function(y), diff_eq_grid, start_value, α=1)

using Plots
plot(map(y -> y[1], solution), map(y -> y[2], solution))
png("plots/task_7_phase_trajectory.png")