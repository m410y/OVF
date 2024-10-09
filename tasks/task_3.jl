include("../src/OVF.jl")
using .OVF
using Plots

f(x) = 1 / (1 + x^2)
a, b = -1, 1
value_true = atan(b) - atan(a)
methods = [
    ("left", OVF.left_integral),
    ("right", OVF.right_integral),
    ("trapezoid", OVF.trapezoid_integral),
    ("middle", OVF.middle_integral),
    ("simpson", OVF.simpson_integral),
]
steps = (b - a) * 2.0 .^ (0:-1:-8)
plot(xaxis = :log, yaxis = :log)
for (method_n, (name, method)) in enumerate(methods)
    values = map(h -> method(f, a:h:b), steps)
    errors = max.(abs.(values .- value_true), eps(value_true))
    plot!(steps, errors, label = name)
    shift, slope = [one.(steps) log.(steps)] \ log.(errors)
    coeff = exp(shift) * sign(values[end] - value_true)
    println("$(name) ~ $(coeff) ⋅ h^ $(slope)")
end
png("plots/task_3_error_order.png")

using SpecialFunctions
plot(x -> OVF.erf(x) - erf(x))
png("plots/task_3_erf_difference.png")
