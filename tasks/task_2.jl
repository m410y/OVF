include("../src/OVF.jl")
using .OVF

U_0 = 1.0
a = 1.0

sym_f(E) = tan(√(2(U_0 + E)) * a) - √(-E / (U_0 + E))
sym_f_diff(E) = cos(√(2(U_0 + E)) * a)^-2 * a / 2√(2(U_0 + E)) + U_0 / 2(U_0 + E) / √(-E * (U_0 + E))

println("sym energy")
println(OVF.dichotomy_root(sym_f, -U_0 + eps(U_0), 0.0, verbose=true))
println(OVF.iteration_root(sym_f, -U_0/2, verbose=true, λ=0.1))
println(OVF.newton_root(sym_f, sym_f_diff, -U_0/2, verbose=true))

println(OVF.dichotomy_root(OVF.j1, 1.0, 5.0, error=1e-15))

# asym_f(E) = cot(√(2(U_0 + E)) * a) - √(-E / (U_0 + E))
# asym_f_diff(E) = -sin(√(2(U_0 + E)) * a)^-2 * a / 2√(2(U_0 + E)) + U_0 / 2(U_0 + E) / √(-E * (U_0 + E))

# println("asym energy")
# println(OVF.dichotomy_root(asym_f, -U_0 + eps(U_0), 0.0, verbose=true))
# println(OVF.iteration_root(asym_f, -U_0/2, verbose=true, λ=0.2))
# println(OVF.newton_root(asym_f, sym_f_diff, -U_0/2, verbose=true))