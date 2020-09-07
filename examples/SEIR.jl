using BenchmarkTools
include("../experiments_bounds.jl")

varnames = ["S", "E", "I", "N", "c", "a", "b", "d", "nu"]
R, (S, E, I, N, c, a, b, d, nu) = PolynomialRing(QQ, varnames)

ode = ODE(Dict(
    S => -b * S * I // N,
    E => b * S * I // N - nu * E,
    I => nu * E - a * I,
    N => R(0),
    c => R(0)
), [])

# compiles the function and runs a number of times, outputs the mean runtime
@btime bound_number_experiments(ode, [c * I + d * E, c, N])

print(bound_number_experiments(ode, [c * I + d * E, c, N]), "\n")
