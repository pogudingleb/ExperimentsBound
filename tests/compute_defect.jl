@testset "Computing identifiability defect" begin

    test_cases = []
    P = fmpq_mpoly
    DType = Union{P, Generic.Frac{P}}

    R, (x, a) = PolynomialRing(QQ, ["x", "a"])
    push!(
        test_cases,
        Dict(
             "ODE" => ODE{P}(Dict{P, DType}(x => x + a), Array{P, 1}()),
            "outputs" => [x],
            "answer" => 0
        )
    )

    R, (x, a) = PolynomialRing(QQ, ["x", "a"])
    push!(
        test_cases,
        Dict(
            "ODE" => ODE{P}(Dict{P, DType}(x => x), Array{P, 1}()),
            "outputs" => [x^2],
            "answer" => 1
        )
    )

    R, (x, y, a, b) = PolynomialRing(QQ, ["x", "y", "a", "b"])
    push!(
        test_cases,
        Dict(
            "ODE" => ODE{P}(Dict{P, DType}(x => a * y, y => b * x), Array{P, 1}()),
            "outputs" => [x],
            "answer" => 1
        )
    )

    R, (x, y, a, b, c, d, u) = PolynomialRing(QQ, ["x", "y", "a", "b", "c", "d", "u"])
    push!(
        test_cases,
        Dict(
             "ODE" => ODE{P}(Dict{P, DType}(
                 x => a * x - b * x * y + u,
                 y => -c * y + d * x * y
             ), [u]),
             "outputs" => [x],
             "answer" => 1
        )
    )

    R, (M, P0, P1, P2, PN, vs, vd, K1, K2, K3, K4, Kd, KI, V1, V2, V3, V4, vm, Km, ks, k1, k2) = PolynomialRing(QQ, ["M", "P0", "P1", "P2", "PN", "vs", "vd", "K1", "K2", "K3", "K4", "Kd", "KI", "V1", "V2", "V3", "V4", "vm", "Km", "ks", "k1", "k2"])
    push!(
        test_cases,
        Dict(
             "ODE" => ODE{P}(Dict{P, DType}(
                 M => vs * KI^4 // (KI^4 + PN^4) - vm * M // (Km + M),
                 P0 => ks * M  - V1 * P0 // (K1 + P0) + V2 * P1 // (K2 + P1),
                 P1 => V1 * P0 // (K1 + P0) + V4 * P2 // (K2 + P2) - P1 * (V2 // (K2 + P1) + V3 // (K3 + P1)),
                 P2 => V3 * P1 // (K3 + P1) - P2 * (V4 // (K4 + P2) + k1 + vd // (Kd + P2)) + k2 * PN,
                 PN => k1 * P2 - k2 * PN
             ), Array{P, 1}()),
             "outputs" => [PN],
             "answer" => 1
        )
    )

    R, (x, a1, a2, a3, a4) = PolynomialRing(QQ, ["x", "a1", "a2", "a3", "a4"])
    push!(
        test_cases,
        Dict(
             "ODE" => ODE{P}(Dict{P, DType}(
                 x => x + a1 + a2 + a3 + a4
             ), Array{P, 1}()),
             "outputs" => [x^5],
             "answer" => 3
        )
    )


    for case in test_cases
        ode, outputs = case["ODE"], case["outputs"]
        @time defect = compute_defect(ode, outputs)
        @test defect == case["answer"]
    end

end
