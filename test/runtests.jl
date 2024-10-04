test_names = [
    "root_find_test",
    "integration_test",
    "specials_test",
    "derivation_test",
    "inter_poly_test",
    "linear_algebra_test",
    "diff_equation_test",
]

for test_name in test_names
    println("Running $(test_name):")
    include("$(test_name).jl")
end
