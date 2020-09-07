using BenchmarkTools
include("../experiments_bounds.jl")

R, (x1, x2, mu1, mu2) = PolynomialRing(QQ, ["x1", "x2", "mu1", "mu2"])

ode = ODE(Dict(
    x1 => R(0),
    x2 => x1 * x2 + mu1 * x1 + mu2
), [])

# compiles the function and runs a number of times, outputs the mean runtime
@btime bound_number_experiments(ode, [x2])

print(bound_number_experiments(ode, [x2]), "\n")
