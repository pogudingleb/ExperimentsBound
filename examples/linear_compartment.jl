function linear_compartment_model(graph, sinks)
    """
    Input: 
        - graph - graph of the model network represented via adjacency lists
        - sinks - the indices of nodes having a sink
    Output: the corresponding ODE object where each parameter a_ij is replaced
            with a_ij + b_ij * x_0, where x_0 is a constant input encoded as a constant output
    """
    n = length(graph)
    x_vars_names = ["x$i" for i in 0:n]
    edges_vars_names = Array{String, 1}()
    for i in 1:n
        for j in graph[i]
            push!(edges_vars_names, "a_$(i)_$(j)")
            push!(edges_vars_names, "b_$(i)_$(j)")
        end
    end
    for s in sinks
        push!(edges_vars_names, "a_$(s)_0")
        push!(edges_vars_names, "b_$(s)_0")
    end
    R, vars = PolynomialRing(QQ, vcat(x_vars_names, edges_vars_names))
    x_vars = vars[2:(n + 1)]
    x0 = vars[1]
    equations = Dict{fmpq_mpoly, Union{fmpq_mpoly, Generic.Frac{fmpq_mpoly}}}(x => R(0) for x in x_vars)
    equations[x0] = R(0)
    for i in 1:n
        for j in graph[i]
            rate = str_to_var("a_$(i)_$(j)", R) + str_to_var("b_$(i)_$(j)", R) * x0

            if i != j
                equations[x_vars[j]] += x_vars[i] * rate
                equations[x_vars[i]] -= x_vars[i] * rate
            else
                equations[x_vars[i]] -= x_vars[i] * rate
            end
        end
        if i in sinks
            rate = str_to_var("a_$(i)_0", R) + str_to_var("b_$(i)_0", R) * x0
            equations[x_vars[i]] += -x_vars[i] * rate
        end
    end
    return ODE{fmpq_mpoly}(equations)
end

#------------------------------------------------------------------------------

function bicycle(n)
    """
    Generates a bidirected cycle of length n
    """
    graph = []
    for i in 1:n
        prev = (i == 1) ? n : (i - 1)
        next = (i == n) ? 1 : i + 1
        push!(graph, [prev, next])
    end
    return graph
end

#------------------------------------------------------------------------------

function cycle(n)
    """
    Single directed cycle
    """
    graph = [[(i == n) ? 1 : (i + 1)] for i in 1:n]
    return graph
end

#------------------------------------------------------------------------------

function catenary(n)
    """
    Bidirected chain from 1 to n
    """
    graph = [[] for i in 1:n]
    for i in 1:n
        if i != 1
            push!(graph[i], i - 1)
        end
        if i != n
            push!(graph[i], i + 1)
        end
    end
    return graph
end

#------------------------------------------------------------------------------

function mammilary(n)
    """
    Bidirected 'star' with center at 1 and rays to 2, ..., n
    """
    graph = []
    push!(graph, [i for i in 2:n])
    for i in 2:n
        push!(graph, [1])
    end
    return graph
end

#------------------------------------------------------------------------------
