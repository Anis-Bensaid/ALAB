function get_path_cost_matrix(paths::Array{Any,1})::Array{Float64, 2}
    C = zeros(num_locomotives, K)
    for i=1:num_locomotives
        for k=1:K
            try
                C[i,k] = paths[i][k][end].f
            catch BoundsError
                C[i,k] = 1e6
            end
        end
    end
    return C
end


function calculate_interaction_penalty(
    path1::Array{Any,1},
    path2::Array{Any,1},
    α::Float64,
    β::Float64
)::Float64
    # Count the number of times the paths cross
    nodes1 = [node.name for node in path1]
    nodes2 = [node.name for node in path2]
    node_crossover_count = max(sum([node in nodes1 for node in nodes2]), sum([node in nodes2 for node in nodes1]))

    # Count the number of edges shared by the paths
    edges1 = [(nodes1[i], nodes1[i+1]) for i=1:length(nodes1)-1]
    edges2 = [(nodes2[i], nodes2[i+1]) for i=1:length(nodes2)-1]
    edge_overlap_count = max(sum([edge in edges1 for edge in edges2]), sum([edge in edges2 for edge in edges1]))

    return α*(node_crossover_count-edge_overlap_count) + β*edge_overlap_count
end


function get_interaction_penalty_array(
    paths::Array{Any,1},
    α::Float64, # Node penalty
    β::Float64 # Edge penalty;
)::Array{Float64, 4}
    P = zeros(num_locomotives, num_locomotives, K, K)
    for i=1:num_locomotives
        for j=1:num_locomotives
            for k=1:K
                for l=1:K
                    try
                        if i != j
                            P[i,j,k,l] = calculate_interaction_penalty(paths[i][k], paths[j][l], α, β)
                        end
                    catch BoundsError
                        P[i,j,k,l] = 1e6
                    end
                end
            end
        end
    end
    return P
end


function optimise_paths(
    paths::Array{Any,1},
    α::Float64 = 1.0, # Node crossover penalty
    β::Float64 = 1.0, # Edge crossover penalty
    γ::Float64 = 1e-3, # Path distance penalty
    solver_time_limit::Int64 = 10
)::Array{Tuple{Int64, Int64}}
    C = get_path_cost_matrix(paths)
    P = get_interaction_penalty_array(paths, α, β)
    num_locomotives, num_paths = size(C)

    model = Model(solver=GurobiSolver(OutputFlag=0, TimeLimit=solver_time_limit, GUROBI_ENV))

    @variable(model, x[1:num_locomotives, 1:num_paths], Bin)
    @variable(model, z[1:num_locomotives, 1:num_locomotives, 1:num_paths, 1:num_paths], Bin)

    @constraint(model, [i=1:num_locomotives], sum(x[i,:]) == 1)
    @constraint(model, [i=1:num_locomotives, j=1:num_locomotives, k=1:num_paths, l=1:num_paths], z[i,j,k,l] <= x[i,k])
    @constraint(model, [i=1:num_locomotives, j=1:num_locomotives, k=1:num_paths, l=1:num_paths], z[i,j,k,l] <= x[j,l])
    @constraint(model, [i=1:num_locomotives, j=1:num_locomotives, k=1:num_paths, l=1:num_paths],
        z[i,j,k,l] >= x[i,k] + x[j,l] - 1)

    @objective(model, Min, sum(P.*z) + γ*sum(C.*x))

    solve(model)
    return Tuple.(findall(getvalue(x) .== 1))
end
