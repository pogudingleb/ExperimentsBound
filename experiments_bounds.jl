using IterTools
using LinearAlgebra
using Logging
using Primes

include("ODE.jl")

#------------------------------------------------------------------------------

function eval_frac(f, point)
    if applicable(numerator, f)
        return divexact(evaluate(numerator(f), point), evaluate(denominator(f), point))
    end
    return evaluate(f, point)
end

#------------------------------------------------------------------------------

function differentiate_solution(ode::ODE, params, ic, inputs, prec::Int)
    """
    Input: the same as for power_series_solutions
    Output: a tuple consisting of the power series solution and 
            a dictionary of the form (u, v) => power series, where u is a state variable
            v is a state of parameter, and the power series is the partial derivative of
            the function u w.r.t. v evaluated at the solution
    """
    @debug "Computing the power series solution of the system"
    ps_sol = power_series_solution(ode, params, ic, inputs, prec)
    ps_ring = parent(Iterators.peel(values(ps_sol))[1])
    for p in ode.parameters
        ps_sol[p] = ps_ring(params[p])
    end
    sol_vec = [ps_sol[v] for v in gens(ode.poly_ring)]

    @debug "Building the variational system at the solution"
    # Y' = AY + B
    vars = vcat(ode.x_vars, ode.parameters)
    SA = MatrixSpace(ps_ring, length(ode.x_vars), length(ode.x_vars))
    A = SA([
        eval_frac(derivative(ode.equations[vars[i]], vars[j]), sol_vec)
        for i in 1:length(ode.x_vars), j in 1:length(ode.x_vars)
    ])
    SB = MatrixSpace(ps_ring, length(ode.x_vars), length(vars))
    B = zero(SB)
    for i in 1:length(ode.x_vars)
        for j in (length(ode.x_vars) + 1):length(vars)
            B[i, j] = eval_frac(derivative(ode.equations[vars[i]], vars[j]), sol_vec)
        end
    end
    # TODO: make use of one() function (problems modulo prime)
    initial_condition = zero(MatrixSpace(base_ring(ode.poly_ring), length(ode.x_vars), length(vars)))
    for i in 1:length(ode.x_vars)
        initial_condition[i, i] = 1
    end
    
    @debug "Solving the variational system and forming the output"
    sol_var_system = ps_matrix_linear_de(A, B, initial_condition, prec)
    return (
        ps_sol, 
        Dict(
            (vars[i], vars[j]) => sol_var_system[i, j]
            for i in 1:length(ode.x_vars), j in 1:length(vars)
        )
    )
end

#------------------------------------------------------------------------------

function differentiate_output(ode::ODE, outputs, params, ic, inputs, prec::Int)
    """
    Similar to differentiate_solution but computes partial derivatives of a prescribed outputs
    returns a list of dictionaries var => dy / dvar
    """
    @debug "Computing partial derivatives of the solution"
    ps_sol, sol_diff = differentiate_solution(ode, params, ic, inputs, prec)
    ps_ring = parent(Iterators.peel(values(ps_sol))[1])
    for p in ode.parameters
        ps_sol[p] = ps_ring(params[p])
    end
    sol_vec = [ps_sol[v] for v in gens(ode.poly_ring)]

    @debug "Evaluating the partial derivatives of the outputs"
    result = []
    for y in outputs
        push!(result, Dict())
        for x in ode.x_vars
            result[end][x] = sum([eval_frac(derivative(y, xx), sol_vec) * sol_diff[(xx, x)] for xx in ode.x_vars])
        end
        for p in ode.parameters
            result[end][p] = sum([eval_frac(derivative(y, xx), sol_vec) * sol_diff[(xx, p)] for xx in ode.x_vars])
            result[end][p] += eval_frac(derivative(y, p), sol_vec)
        end
    end

    return result 
end

#------------------------------------------------------------------------------

function get_degree_and_coeffsize(f)
    """
    for f being rational function or polynomial over QQ returns a tuple
    (degree, max_coef_size)
    """
    if applicable(numerator, f)
        num_deg, num_coef = get_degree_and_coeffsize(numerator(f))
        den_deg, den_coef = get_degree_and_coeffsize(denominator(f))
        return (max(num_deg, den_deg), max(num_coef, den_coef))
    end
    if length(f) == 0
        return (0, 1)
    end
    max_coef = 1
    for c in coeffs(f)
        max_coef = max(max_coef, 2 * height_bits(c))
    end
    return (total_degree(f), max_coef)
end

#------------------------------------------------------------------------------

function compute_defect(ode::ODE, outputs, p::Float64 = 0.99)
    """
    Computed the identifiability defect (Definition 2.6) of an ode system with respect
    to the given outputs with probability at least p
    """

    # Computing the prime using Theorem 1.1 from https://doi.org/10.1006/jsco.2002.0532
    @debug "Computing the prime number"
    d, h = 1, 1
    for f in vcat(collect(values(ode.equations)), outputs)
        df, hf = get_degree_and_coeffsize(f)
        d = max(d, df)
        h = max(h, hf)
    end
    mu = ceil(1 / (1 - sqrt(p)))
    n = length(ode.x_vars)
    m = length(outputs)
    r = length(ode.u_vars)
    ell = length(ode.parameters)
    D = 4 * (n + ell)^2 * (n + m) * d
    Dprime = D * (2 * log(n + ell + r + 1) + log(mu * D)) + 4 * (n + ell)^2 * ((n + m) * h + log(2 * n * D))
    prime = Primes.nextprime(Int(ceil(2 * mu * Dprime)))
    @debug "The prime is $prime"
    F = GF(prime)
 
    @debug "Reducing the system modulo prime"
    ode_red = reduce_ode_mod_p(ode, prime)
    outputs_red = [_reduce_poly_mod_p(poly, prime) for poly in outputs]
    prec = length(ode.x_vars) + length(ode.parameters)
    params_vals = Dict(p => F(rand(1:prime)) for p in ode_red.parameters)
    ic = Dict(x => F(rand(1:prime)) for x in ode_red.x_vars)
    inputs = Dict(u => [F(rand(1:prime)) for i in 1:prec] for u in ode_red.u_vars)

    @debug "Computing the output derivatives"
    output_derivatives = differentiate_output(ode_red, outputs_red, params_vals, ic, inputs, prec)

    @debug "Building the matrices"
    Jac_x_param = zero(MatrixSpace(F, length(outputs_red) * prec, prec))
    xs_params = vcat(ode_red.x_vars, ode_red.parameters)
    for i in 1:length(outputs_red)
        for j in 1:prec
            for k in 1:length(xs_params)
                Jac_x_param[(i - 1) * prec + j, k] = coeff(output_derivatives[i][xs_params[k]], j - 1)
            end
        end
    end
    Jac_x = Jac_x_param[:, 1:length(ode_red.x_vars)]

    @debug "Computing the ranks"
    return LinearAlgebra.rank(Jac_x) + ell - LinearAlgebra.rank(Jac_x_param)
end

#------------------------------------------------------------------------------

function bound_number_experiments(ode::ODE, outputs, p::Float64=0.99)
    """
    Input:
        - ode, an ODE object representing an ODE system over Q
        - outputs, rational function being the outputs of the system
    Output: a named tuple with field loc and glob giving upper bounds to the
    NumExpLoc and NumExpGlob, respectively (for definitions, see the paper)
    """
    defect_prev = length(ode.parameters)
    for i in 1:(length(ode.parameters) + 1)
        ode_replicated, outputs_replicated = generate_replica(ode, outputs, i)
        defect_cur = compute_defect(ode_replicated, outputs_replicated, 1. - (1. - p) / length(ode.parameters))
        if defect_prev == defect_cur
            return (loc = i - 1, glob = i)
        end
        defect_prev = defect_cur
    end
    throw(Core.ErrorException("Too many steps for defect to stabilize. Something is wrong here"))
end

#------------------------------------------------------------------------------
