# Gleb: the code taken from the ode elimination projects, tests are there

using Oscar

include("power_series_utils.jl")

struct ODE
    poly_ring
    x_vars
    u_vars
    parameters
    equations
    
    function ODE(eqs, inputs)
        #Initialize ODE
        #equations is a dictionary x_i => f_i(x, u, params)

        num, den = unpack_fraction(collect(values(eqs))[1])
        poly_ring = parent(num)
        x_vars = collect(keys(eqs))
        u_vars = inputs
        parameters = filter(v -> (!(v in x_vars) && !(v in u_vars)), gens(poly_ring))
        new(poly_ring, x_vars, u_vars, parameters, eqs)
    end
end

#------------------------------------------------------------------------------

function power_series_solution(ode::ODE, param_values, initial_conditions, input_values, prec)
    """
    Input:
        - ode, an ode to solve
        - param_values, initial_conditions - parameter values and initial conditions to plug in
          both are dictionaries variable => value
        - input_values - power series for the inpiuts presented as a dictionary
          variable => list of coefficients
        - prec, the precision
    Output: computes a power series solution with precision prec presented as a dictionary
            variable => corresponding coordiante of the solution
    """
    new_varnames = map(string, vcat(ode.x_vars, map(v -> "$(v)_dot", ode.x_vars), ode.u_vars))

    new_ring, new_vars = PolynomialRing(base_ring(ode.poly_ring), new_varnames)
    equations = Array{RingElem, 1}()
    evaluation = Dict(k => new_ring(v) for (k, v) in param_values)
    for v in vcat(ode.x_vars, ode.u_vars)
        evaluation[v] = str_to_var(string(v), new_ring)
    end
    for (v, eq) in ode.equations
        num, den = map(p -> eval_at_dict(p, evaluation), unpack_fraction(eq))
        push!(equations, den * str_to_var("$(v)_dot", new_ring) - num)
    end
    new_inputs = Dict(str_to_var(string(k), new_ring) => v for (k, v) in input_values)
    new_ic = Dict(str_to_var(string(k), new_ring) => v for (k, v) in initial_conditions)
    result = ps_ode_solution(equations, new_ic, new_inputs, prec)
    return Dict(v => result[str_to_var(string(v), new_ring)] for v in vcat(ode.x_vars, ode.u_vars))
end

#------------------------------------------------------------------------------

function _reduce_poly_mod_p(poly, p)
    """
    Reduces a polynomial over Q modulo p
    """
    den = denominator(poly)
    num = change_base_ring(ZZ, den * poly)
    if GF(p)(den) == 0
        throw(Base.ArgumentError("Prime $p divides the denominator of $poly"))
    end
    return change_base_ring(GF(p), num) * (1 // GF(p)(den))
end

#--------------------------------------

function reduce_ode_mod_p(ode::ODE, p)
    """
    Input: ode is an ODE over QQ, p is a prime number
    Output: the reduction mod p, throws an exception if p divides one of the denominators
    """
    new_ring, new_vars = PolynomialRing(GF(p), map(string, gens(ode.poly_ring)))
    new_inputs = map(u -> str_to_var(string(u), new_ring), ode.u_vars)
    new_eqs = Dict()
    for (v, f) in ode.equations
        new_v = str_to_var(string(v), new_ring)
        if applicable(numerator, f)
            # if f is a rational function
            num, den = map(poly -> _reduce_poly_mod_p(poly, p), [numerator(f), denominator(f)])
            if den == 0
                throw(Base.ArgumentError("Prime $p divides the denominator of $poly"))
            end
            new_eqs[new_v] = num // den
        else
            new_eqs[new_v] = _reduce_poly_mod_p(f, p)
        end
    end
    return ODE(new_eqs, new_inputs)
end

#------------------------------------------------------------------------------

function generate_replica(ode::ODE, outputs, r::Int)
    """
    Returns (ode_r, outputs_r), and r-fold replica of the original pair (ode, outputs).
    States, outputs, and inputs are replicated, parameters are not
    """
    new_varnames = map(string, ode.parameters)
    for v in vcat(ode.x_vars, ode.u_vars)
        append!(new_varnames, ["$(v)_r$i" for i in 1:r])
    end
    new_ring, new_vars = PolynomialRing(base_ring(ode.poly_ring), new_varnames)
    new_eqs = Dict()
    new_outputs = []
    new_us = []
    for i in 1:r
        eval = merge(
            Dict(p => str_to_var(string(p), new_ring) for p in ode.parameters),
            Dict(v => str_to_var("$(v)_r$i", new_ring) for v in vcat(ode.x_vars, ode.u_vars))
        )
        eval_vec = [eval[v] for v in gens(ode.poly_ring)]
        append!(new_outputs, map(p -> evaluate(p, eval_vec), outputs))
        new_eqs = merge(
            new_eqs, 
            Dict(evaluate(x, eval_vec) => evaluate(f, eval_vec) for (x, f) in ode.equations)
        )
        append!(new_us, [str_to_var("$(u)_r$i", new_ring) for u in ode.u_vars])
    end
    return (ODE(new_eqs, new_us), new_outputs)
end

#------------------------------------------------------------------------------
