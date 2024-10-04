include("../src/OVF.jl")
using .OVF

atan_integral(method; step = 1e-1) = begin
    int = method(x -> 1 / (1 + x^2), -1:step:1)
    println("$(method): $(int), $(int - pi/2)")
end

foreach(
    atan_integral,
    [
        OVF.left_integral,
        OVF.right_integral,
        OVF.trapezoid_integral,
        OVF.middle_integral,
        OVF.simpson_integral,
    ],
)

using Plots
using SpecialFunctions
plot(x -> OVF.erf(x) - erf(x))
png("plots/task_3_erf_difference.png")
