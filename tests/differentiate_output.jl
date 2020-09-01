@testset "Partial derivatives of a solution w.r.t. to initial conditions and parameters" begin

    test_cases = []

    R, (x, a) = PolynomialRing(QQ, ["x", "a"])
    push!(
        test_cases,
        Dict(
          "ODE" => ODE(Dict(x => x + a), []),
          "ic" => Dict(x => rand(1:10)),
          "param_vals" => Dict(a => rand(1:10)),
          "inputs" => Dict(),
          "prec" => 20,
          "outputs" => [x^2]
        )
    )

    R, (x, a) = PolynomialRing(QQ, ["x", "a"])
    push!(
        test_cases,
        Dict(
          "ODE" => ODE(Dict(x => x^2 + a), []),
          "ic" => Dict(x => rand(1:10)),
          "param_vals" => Dict(a => rand(1:10)),
          "inputs" => Dict(),
          "prec" => 20,
          "outputs" => [x + a^2, x^3]
        )
    )

    R, (x, y, a, b) = PolynomialRing(QQ, ["x", "y", "a", "b"])
    push!(
        test_cases,
        Dict(
          "ODE" => ODE(Dict(x => x^2 + 2 * x * y - 3 * a * y, y =>  x^2 + a * b - b^2 + 4 * b * x), []),
          "ic" => Dict(x => rand(1:10), y => rand(1:10)),
          "param_vals" => Dict(a => rand(1:10), b => rand(1:10)),
          "inputs" => Dict(),
          "prec" => 8,
          "outputs" => [a * x, b * y^2 - y]
        )
    )

    R, (x, a, u) = PolynomialRing(QQ, ["x", "a", "u"])
    push!(
        test_cases,
        Dict(
          "ODE" => ODE(Dict(x => u + a), [u]),
          "ic" => Dict(x => rand(1:10)),
          "param_vals" => Dict(a => rand(1:10)),
          "inputs" => Dict(u => [rand(-3:3) for i in 1:20]),
          "prec" => 20,
          "outputs" => [x]
        )
    )

    varnames = vcat(
        ["x_$i" for i in 1:3],
        ["p_$i" for i in 1:3],
        ["u_$i" for i in 1:2],
    )
    R, vars = PolynomialRing(GF(2^31 - 1), varnames)
    push!(
        test_cases,
        Dict(
          "ODE" => ODE(Dict(vars[i] => rand_poly(1, vars) for i in 1:3), vars[7:8]),
          "ic" => Dict(vars[i] => rand(1:50) for i in 1:3),
          "param_vals" => Dict(vars[i + 3] => rand(1:50) for i in 1:3),
          "inputs" => Dict(u => [rand(-30:30) for i in 1:6] for u in vars[7:end]),
          "prec" => 6,
          "outputs" => [rand_poly(2, vars) for i in 1:3]
        )
    )

    varnames = vcat(
        ["x_$i" for i in 1:3],
        ["p_$i" for i in 1:3],
        ["u_$i" for i in 1:2],
    )
    R, vars = PolynomialRing(GF(2^31 - 1), varnames)
    push!(
        test_cases,
        Dict(
          "ODE" => ODE(Dict(vars[i] => rand_poly(2, vars) for i in 1:3), vars[7:8]),
          "ic" => Dict(vars[i] => rand(1:50) for i in 1:3),
          "param_vals" => Dict(vars[i + 3] => rand(1:50) for i in 1:3),
          "inputs" => Dict(u => [rand(-30:30) for i in 1:6] for u in vars[7:end]),
          "prec" => 6,
          "outputs" => [rand_poly(2, vars) for i in 1:3]
        )
    )
 
    varnames = vcat(
        ["x_$i" for i in 1:2],
        ["p_$i" for i in 1:2],
        "u",
    )
    R, vars = PolynomialRing(GF(2^31 - 1), varnames)
    push!(
        test_cases,
        Dict(
          "ODE" => ODE(Dict(vars[i] => rand_poly(1, vars) // (vars[1] + vars[3]) for i in 1:2), [vars[end]]),
          "ic" => Dict(vars[i] => rand(1:50) for i in 1:2),
          "param_vals" => Dict(vars[i + 2] => rand(1:50) for i in 1:2),
          "inputs" => Dict(vars[end] => [rand(-30:30) for i in 1:4]),
          "prec" => 4,
          "outputs" => [rand_poly(1, vars) for i in 1:2]
        )
    )
 
    for case in test_cases
        ode, prec = case["ODE"], case["prec"]
        outputs = case["outputs"]
        @time sol1 = differentiate_output(ode, outputs, case["param_vals"], case["ic"], case["inputs"], prec)
        sol2 = diff_sol_Lie_derivatives(ode, outputs, case["param_vals"], case["ic"], case["inputs"], prec)
        for i in 1:length(outputs)
            for v in vcat(ode.x_vars, ode.parameters)
                @test sol2[i][v] == [base_ring(ode.poly_ring)(coeff(sol1[i][v], j) * factorial(j)) for j in 0:(prec - 1)]
            end
        end
    end

end
