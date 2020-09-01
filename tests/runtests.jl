using Test
using TestSetExtensions

using Oscar

include("../experiments_bounds.jl")

#------- Auxiliary functions --------------------------------------------------

function diff_sol_Lie_derivatives(ode::ODE, outputs, params, ic, inputs, prec::Int)
    # creating a new ring with variables for the derivatives of u
    new_varnames = map(string, gens(ode.poly_ring))
    if length(ode.u_vars) > 0
        append!(new_varnames, ["$(u)_$i" for u in ode.u_vars for i in 1:(prec - 1)])
    end
    new_ring, vars = PolynomialRing(base_ring(ode.poly_ring), new_varnames)
    
    # mapping everything to the new ring
    eval_point = [str_to_var(string(v), new_ring) for v in gens(ode.poly_ring)]
    new_eqs = Dict()
    for (x, f) in ode.equations
        if applicable(f, numerator)
            new_eqs[str_to_var(string(x), new_ring)] = evaluate(numerator(f), eval_point) // evaluate(denominator(f), eval_point)
        else
            new_eqs[str_to_var(string(x), new_ring)] = evaluate(f, eval_point)
        end
    end
    outputs_new = []
    for g in outputs
        if applicable(g, numerator)
            push!(outputs_new, evaluate(numerator(g), eval_point) // evaluate(denominator(g), eval_point))
        else
            push!(outputs_new, evaluate(g, eval_point))
        end
    end
    outputs = outputs_new
    params, ic = map(d -> Dict(str_to_var(string(k), new_ring) => v for (k, v) in d), [params, ic])

    # computing Lie derivatives
    derivation = copy(new_eqs)
    for u in ode.u_vars
        derivation[str_to_var(string(u), new_ring)] = str_to_var("$(u)_1", new_ring)
        for i in 1:(prec - 2)
            derivation[str_to_var("$(u)_$i", new_ring)] = str_to_var("$(u)_$(i + 1)", new_ring)
        end
    end
    Lie_derivatives = []
    for g in outputs
        push!(Lie_derivatives, [g // new_ring(1)])
        for i in 1:prec
            push!(
                Lie_derivatives[end],
                sum([derivative(Lie_derivatives[end][end], v) * get(derivation, v, 0) for v in gens(new_ring)])
            )
        end
    end

    # producing the result
    eval_dict = merge(params, ic)
    for u in ode.u_vars
        eval_dict[str_to_var(string(u), new_ring)] = inputs[u][1]
        for i in 2:prec
            eval_dict[str_to_var("$(u)_$(i - 1)", new_ring)] = inputs[u][i] * factorial(i - 1)
        end
    end
    eval_vec = [eval_dict[v] for v in gens(new_ring)]

    result = []
    for i in 1:length(outputs)
        push!(result, Dict())
        for v in vcat(ode.x_vars, ode.parameters)
            result[end][v] = []
            for j in 1:prec
                push!(
                    result[end][v], 
                    eval_frac(derivative(Lie_derivatives[i][j], str_to_var("$v", new_ring)), eval_vec)
                )
            end
        end
    end

    return result
end

#------------------------------------------------------------------------------

function rand_poly(deg, vars)
    if deg == 0
        return parent(vars[1])(1)
    end
    result = 0
    indices = collect(1:length(vars))
    monomials = []
    for d in 0:deg
        for subs in IterTools.subsets(indices, d)
            push!(monomials, subs)
        end
    end

    for subs in monomials
        monom = rand(-50:50)
        for v_ind in subs
            monom *= vars[v_ind]
        end
        result += monom
    end

    return result
end

#------------------------------------------------------------------------------

@info "Testing started"

@testset "All the tests" begin
    @includetests ARGS
end

