using BenchmarkTools
include("../experiments_bounds.jl")

ode = @ODEmodel(
    S'(t) = -b * S(t) * I(t) / N(t),
    E'(t) = b * S(t) * I(t) / N(t) - nu * E(t),
    I'(t) = nu * E(t) - a * I(t),
    N'(t) = 0,
    c'(t) = 0,
    [d]
)

# compiles the function and runs a number of times, outputs the mean runtime
@btime bound_number_experiments(ode, [c * I + d * E, c, N])

print(bound_number_experiments(ode, [c * I + d * E, c, N]), "\n")
