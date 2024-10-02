include("../src/OVF.jl")
using .OVF

diff_eq_grid = range(
    start=-pi/2,
    stop=pi/2,
    length=10
)
A_left, B_left, C_left = 1., -2., 1.
A_right, B_right, C_right = 3., 1., -1.

h = step(diff_eq_grid)
one_line =  fill(1, length(diff_eq_grid) - 2)
solution = OVF.thomas_algorithm!(
    [one_line/h^2; -B_right/h],
    [A_left - B_left/h; -2one_line/h^2; A_right + B_right/h],
    [B_left/h; one_line/h^2],
    [C_left; cos.(diff_eq_grid[begin+1:end-1]); C_right]
)

using Plots
a, b = [A_left (B_left - pi/2*A_left); A_right (B_right + pi/2*A_right)] \ [C_left + B_left, C_right - B_right]
plot(x -> a + b * x - cos(x), xlims=(diff_eq_grid[begin], diff_eq_grid[end]))
scatter!(diff_eq_grid, solution)
png("plots/task_9_solution.jl")