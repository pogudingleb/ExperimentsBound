@testset "Partial derivatives of a solution w.r.t. to initial conditions and parameters" begin

    test_cases = []
    DType = Union{fmpq_mpoly, Generic.Frac{fmpq_mpoly}}

    R, (x, a) = PolynomialRing(QQ, ["x", "a"])
    push!(
        test_cases,
        Dict(
          "ODE" => ODE{fmpq_mpoly}(Dict{fmpq_mpoly, DType}(x => x + a), Array{fmpq_mpoly, 1}()),
          "ic" => Dict(x => QQ(rand(1:10))),
          "param_vals" => Dict(a => QQ(rand(1:10))),
          "inputs" => Dict{fmpq_mpoly, Array{fmpq, 1}}(),
          "prec" => 20
        )
    )

    R, (x, a) = PolynomialRing(QQ, ["x", "a"])
    push!(
        test_cases,
        Dict(
          "ODE" => ODE{fmpq_mpoly}(Dict{fmpq_mpoly, DType}(x => x^2 + a),Array{fmpq_mpoly, 1}()),
          "ic" => Dict(x => QQ(rand(1:10))),
          "param_vals" => Dict(a => QQ(rand(1:10))),
          "inputs" => Dict{fmpq_mpoly, Array{fmpq, 1}}(),
          "prec" => 20
        )
    )

    R, (x, y, a, b) = PolynomialRing(QQ, ["x", "y", "a", "b"])
    push!(
        test_cases,
        Dict(
          "ODE" => ODE{fmpq_mpoly}(Dict{fmpq_mpoly, DType}(x => x^2 + 2 * x * y - 3 * a * y, y =>  x^2 + a * b - b^2 + 4 * b * x), Array{fmpq_mpoly, 1}()),
          "ic" => Dict(x => QQ(rand(1:10)), y => QQ(rand(1:10))),
          "param_vals" => Dict(a => QQ(rand(1:10)), b => QQ(rand(1:10))),
          "inputs" => Dict{fmpq_mpoly, Array{fmpq, 1}}(),
          "prec" => 8
        )
    )

    R, (x, a, u) = PolynomialRing(QQ, ["x", "a", "u"])
    push!(
        test_cases,
        Dict(
          "ODE" => ODE{fmpq_mpoly}(Dict{fmpq_mpoly, DType}(x => u + a), [u]),
          "ic" => Dict(x => QQ(rand(1:10))),
          "param_vals" => Dict(a => QQ(rand(1:10))),
          "inputs" => Dict(u => [QQ(rand(-3:3)) for i in 1:20]),
          "prec" => 20
        )
    )

    F = GF(2^31 - 1)
    P = Generic.MPoly{Nemo.gfp_elem}
    DType = Union{P, Generic.Frac{P}}

    varnames = vcat(
        ["x_$i" for i in 1:3],
        ["p_$i" for i in 1:3],
        ["u_$i" for i in 1:2],
    )
    R, vars = PolynomialRing(F, varnames)
    push!(
        test_cases,
        Dict(
          "ODE" => ODE{P}(Dict{P, DType}(vars[i] => rand_poly(1, vars) for i in 1:3), vars[7:8]),
          "ic" => Dict(vars[i] => F(rand(1:50)) for i in 1:3),
          "param_vals" => Dict(vars[i + 3] => F(rand(1:50)) for i in 1:3),
          "inputs" => Dict(u => [F(rand(-30:30)) for i in 1:6] for u in vars[7:end]),
          "prec" => 6
        )
    )

    varnames = vcat(
        ["x_$i" for i in 1:3],
        ["p_$i" for i in 1:3],
        ["u_$i" for i in 1:2],
    )
    R, vars = PolynomialRing(F, varnames)
    push!(
        test_cases,
        Dict(
          "ODE" => ODE{P}(Dict{P, DType}(vars[i] => rand_poly(2, vars) for i in 1:3), vars[7:8]),
          "ic" => Dict(vars[i] => F(rand(1:50)) for i in 1:3),
          "param_vals" => Dict(vars[i + 3] => F(rand(1:50)) for i in 1:3),
          "inputs" => Dict(u => [F(rand(-30:30)) for i in 1:6] for u in vars[7:end]),
          "prec" => 6
        )
    )
 
    varnames = vcat(
        ["x_$i" for i in 1:2],
        ["p_$i" for i in 1:2],
        "u",
    )
    R, vars = PolynomialRing(F, varnames)
    push!(
        test_cases,
        Dict(
          "ODE" => ODE{P}(Dict{P, DType}(vars[i] => rand_poly(1, vars) // (vars[1] + vars[3]) for i in 1:2), [vars[end]]),
          "ic" => Dict(vars[i] => F(rand(1:50)) for i in 1:2),
          "param_vals" => Dict(vars[i + 2] => F(rand(1:50)) for i in 1:2),
          "inputs" => Dict(vars[end] => [F(rand(-30:30)) for i in 1:4]),
          "prec" => 4
        )
    )
 
    for case in test_cases
        ode, prec = case["ODE"], case["prec"]
        @time sol1 = differentiate_solution(ode, case["param_vals"], case["ic"], case["inputs"], prec)[2]
        sol2 = diff_sol_Lie_derivatives(ode, ode.x_vars, case["param_vals"], case["ic"], case["inputs"], prec)
        for i in 1:length(ode.x_vars)
            x = ode.x_vars[i]
            for v in vcat(ode.x_vars, ode.parameters)
                @test sol2[i][v] == [base_ring(ode.poly_ring)(coeff(sol1[(x, v)], i) * factorial(i)) for i in 0:(prec - 1)]
            end
        end
    end

end
