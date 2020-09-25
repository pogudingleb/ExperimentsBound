using BenchmarkTools
include("../experiments_bounds.jl")

ode = @ODEmodel(
    x1'(t) = 0,
    x2'(t) = x1(t) * x2(t) + mu1 * x1(t) + mu2
)

# compiles the function and runs a number of times, outputs the mean runtime
@btime bound_number_experiments(ode, [x2])

print(bound_number_experiments(ode, [x2]), "\n")
