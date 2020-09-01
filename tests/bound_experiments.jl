@testset "Bounding number of experiments" begin

    test_cases = []

    R, (x, a) = PolynomialRing(QQ, ["x", "a"])
    push!(
        test_cases,
        Dict(
            "ODE" => ODE(Dict(x => x + a), []),
            "outputs" => [x],
            "answer" => (loc = 1, glob = 2)
        )
    )

    R, (x, a) = PolynomialRing(QQ, ["x", "a"])
    push!(
        test_cases,
        Dict(
            "ODE" => ODE(Dict(x => x), []),
            "outputs" => [x^2],
            "answer" => (loc = 0, glob = 1)
        )
    )

    R, (x, a) = PolynomialRing(QQ, ["x", "a"])
    push!(
        test_cases,
        Dict(
            "ODE" => ODE(Dict(x => R(0)), []),
            "outputs" => [a^2 * x + a, x],
            "answer" => (loc = 1, glob = 2)
        )
    )

    R, (k1, k2, eA, eB, eC, xA, xB, xC) = PolynomialRing(QQ, ["k1", "k2", "eA", "eB", "eC", "xA", "xB", "xC"])
    push!(
        test_cases,
        Dict(
            "ODE" => ODE(Dict(
                xA => -k1 * xA,
                xB => k1 * xA - k2 * xB,
                xC => k2 * xB,
                eA => R(0),
                eC => R(0)                 
            ), []),
            "outputs" => [xC, eA * xA + eB * xB + eC * xC, eA, eC],
            "answer" => (loc = 1, glob = 2)
        )
    )


    n12, n21, n13, n31, n23, n32 = 15, 13, 12, 14, 17, 20
    R, (a12, a21, a13, a31, a23, a32, x1, x2, x3, u) = PolynomialRing(QQ, ["a12", "a21", "a13", "a31", "a23", "a32", "x1", "x2", "x3", "u"])
    push!(
        test_cases,
        Dict(
            "ODE" => ODE(Dict(
                x1 => a12 * u^n12 * x2 + a13 * u^n13 * x3 - u^n12 * a12 * x1 - a13 * u^n13 * x1,
                x2 => a21 * u^n21 * x1 + a23 * u^n23 * x3 - u^n21 * a21 * x2 - a23 * u^n23 * x2,
                x3 => a31 * u^n31 * x1 + a32 * u^n32 * x2 - u^n31 * a31 * x3 - a32 * u^n32 * x3,
                u => R(0)
            ), []),
            "outputs" => [x1, u],
            "answer" => (loc = 3, glob = 4)
        )
    )

    R, (a, b, x1, x2) = PolynomialRing(QQ, ["a", "b", "x1", "x2"])
    push!(
        test_cases,
        Dict(
            "ODE" => ODE(Dict(
                x1 => x1 * x2 - a * x2 + b,
                x2 => R(0)
            ), []),
            "outputs" => [x1],
            "answer" => (loc = 2, glob = 3)
        )
    )

    for case in test_cases
        ode, outputs = case["ODE"], case["outputs"]
        @time num_exp = bound_number_experiments(ode, outputs)
        @test num_exp == case["answer"]
    end

end
