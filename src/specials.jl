include("integration.jl")

erf(x::Real) = begin
    abs(x) < 6.0 || return sign(x)
    int = 0.0
    if abs(x) < 1.0
        int_range = range(0, abs(x), length=2 + Int(abs(x) ÷ 1e-3))
        int += simpson_integral(t -> exp(-t^2), int_range)
    else
        int += simpson_integral(t -> exp(-t^2), 0:1e-3:1)
        n_bins = Int(abs(x) ÷ 1e-3)
        int_range = range(1.0, abs(x), length=min(9999, n_bins))
        int += simpson_integral(t -> exp(-t^2), int_range)
    end
    return sign(x) * 2/√pi * int
end

j0(x::Real) = 1/pi * simpson_integral(t -> cos(x * sin(t)), range(0, pi, 64))
j1(x::Real) = 1/pi * simpson_integral(t -> cos(t - x * sin(t)), range(0, pi, 64))
jm(x::Real, m::Integer) = 1/pi * simpson_integral(t -> cos(m * t - x * sin(t)), range(0, pi, 64))