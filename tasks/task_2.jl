include("../src/OVF.jl")
using .OVF
using Plots

U_0 = 1.0
a = 1.0

sym_f(E) = tan(√(2(U_0 + E)) * a) - √(-E / (U_0 + E))
sym_f_diff(E) =
    cos(√(2(U_0 + E)) * a)^-2 * a / √(2(U_0 + E)) + U_0 / 2(U_0 + E) / √(-E * (U_0 + E))

println(OVF.dichotomy_root(sym_f, -U_0 + eps(U_0), 0.0, show_iter = true))
println(OVF.iteration_root(sym_f, -U_0 / 2, show_iter = true, λ = 0.1))
println(OVF.newton_root(sym_f, sym_f_diff, -U_0 / 2, show_iter = true))
println(OVF.secant_root(sym_f, -U_0 / 2, show_iter = true, λ = 0.1))
