function linear_compartment_model(graph, sinks, inputs)
    """
    Input: 
        - graph - graph of the model network represented via adjacency lists
        - sinks - the indices of nodes having a sink
        - inputs - list of the nodes having differentially transcendental input
    Output: the corresponding ODE object where each parameter a_ij is replaced
            with a_ij + b_ij * v, where v is a constant input encoded as a constant output
    """
    n = length(graph)
    x_vars_names = ["x$i" for i in 1:n]
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
    u_vars_names = ["u$i" for i in inputs]
    R, vars = PolynomialRing(QQ, vcat(x_vars_names, edges_vars_names, u_vars_names, ["v"]))
    x_vars = vars[1:n]
    v = vars[end]
    equations = Dict(x => R(0) for x in x_vars)
    equations[v] = R(0)
    for i in 1:n
        for j in graph[i]
            rate = str_to_var("a_$(i)_$(j)", R) + str_to_var("b_$(i)_$(j)", R) * v
            if i != j
                equations[x_vars[j]] += x_vars[i] * rate
                equations[x_vars[i]] -= x_vars[i] * rate
            else
                equations[x_vars[i]] -= x_vars[i] * rate
            end
        end
        if i in sinks
            rate = str_to_var("a_$(i)_0", R) + str_to_var("b_$(i)_0", R) * v
            equations[x_vars[i]] += -x_vars[i] * rate
        end
        if i in inputs
            equations[x_vars[i]] += str_to_var("u$i", R)
        end
    end
    return ODE(equations, [str_to_var("u$i", R) for i in inputs])
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
