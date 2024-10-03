include("../src/OVF.jl")
using .OVF

diff_eq_grid = range(
    start=-pi/2,
    stop=pi/2,
    length=10
)
A_left, B_left, C_left = left_borders = 1, -2, 1
A_right, B_right, C_right = right_borders = 3, 1, -1

solution = OVF.poisson_1D(cos, diff_eq_grid, left_borders, right_borders)

using Plots
a, b = [A_left  (B_left - pi/2*A_left);
        A_right (B_right + pi/2*A_right)] \
    [C_left + B_left, C_right - B_right]
plot(x -> a + b * x - cos(x), xlims=(diff_eq_grid[begin], diff_eq_grid[end]))
scatter!(diff_eq_grid, solution)
png("plots/task_9_solution.jl")