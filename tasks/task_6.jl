include("../src/OVF.jl")
using .OVF

diff_eq_grid = range(
    start=0,
    stop=3,
    length=5
)
diff_eq_function = (x, y) -> -y
start_value = 1.0
diff_eq_solution = x -> exp(-x)

euler_solution = OVF.euler_method(diff_eq_function, diff_eq_grid, start_value)
rk2_solution = OVF.rk2_method(diff_eq_function, diff_eq_grid, start_value, α=1)
rk4_solution = OVF.rk4_method(diff_eq_function, diff_eq_grid, start_value)

using Plots
plot(diff_eq_solution, xlims=(diff_eq_grid[begin]-0.1, diff_eq_grid[end]+0.1))
scatter!(diff_eq_grid, euler_solution)
scatter!(diff_eq_grid, rk2_solution)
scatter!(diff_eq_grid, rk4_solution)
png("plots/task_6_solutions.png")