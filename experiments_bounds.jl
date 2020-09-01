using IterTools

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
    ps_sol = power_series_solution(ode, params, ic, inputs, prec)
    ps_ring = parent(Iterators.peel(values(ps_sol))[1])
    for p in ode.parameters
        ps_sol[p] = ps_ring(params[p])
    end
    sol_vec = [ps_sol[v] for v in gens(ode.poly_ring)]

    # building the variational system at the solution
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
    
    # solving the variational system and forming the output
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
    ps_sol, sol_diff = differentiate_solution(ode, params, ic, inputs, prec)
    ps_ring = parent(Iterators.peel(values(ps_sol))[1])
    for p in ode.parameters
        ps_sol[p] = ps_ring(params_vals[p])
    end
    sol_vec = [ps_sol[v] for v in gens(ode.poly_ring)]

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

function compute_defect(ode::ODE, outputs, p::Float64 = 0.99)
    """
    Computed the identifiability defect (Definition 2.6) of an ode system with respect
    to the given outputs with probability at least p
    """

    prime = 2^31 - 1
    F = GF(prime)
 
    ode_red = reduce_ode_mod_p(ode, prime)
    outputs_red = [_reduce_poly_mod_p(poly, prime) for poly in outputs]
    prec = length(ode.x_vars) + length(ode.parameters)
    params_vals = Dict(p => F(rand(1:prime)) for p in ode_red.params)
    ic = Dict(x => F(rand(1:prime)) for x in ode_red.x_vars)
    inputs = Dict(u => [F(rand(1:prime)) for i in 1:prec] for u in ode_red.u_vars)

    output_derivatives = differentiate_output(ode_red, outputs_red, params_vals, ic, inputs, prec)

    # building the matrices
    Jac_x_param = zero(MatrixSpace(F, length(outputs_red) * prec, prec))
    xs_params = vcat(ode_red.x_vars, ode_red.parameters)
    for i in 1:length(outputs_red)
        for j in 1:prec
            for k in 1:length(xs_params)
                Jac_x_param[(i - 1) * prec + j, k] = output_derivatives[i][xs_params[k]][j]
            end
        end
    end
    Jac_x = Jac_x_param

    return LinearAlgebra.rank(Jac_x_params) - LinearAlgebra.rank(Jac_x)
end

#------------------------------------------------------------------------------
