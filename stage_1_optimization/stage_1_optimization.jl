using JuMP, Gurobi, Dates
GUROBI_ENV = Gurobi.Env()

# Make adacency array of network given problem dimensions
function make_adjacency_array(
    num_nodes::Int64,
    num_periods::Int64
)::Array{Int64, 6}
    A = zeros(Int, num_nodes, num_nodes, num_periods, num_periods, 2, 2)

    for t=1:num_periods-1

        if t < num_periods+1
            # The depot pickup node at period t  is connected to all non-depot delivery nodes at period t+1
            A[1, 2:end, t, t+1, 1, 2] .= 1

            # The depot delivery node at period t is connected to all non-depot pickup nodes at period t+1
            A[1, 2:end, t, t+1, 2, 1] .= 1

            # The depot pickup node at period t is connected to the depot pickup node at period t+1
            A[1, 1, t, t+1, 1, 1] = 1
        end

        if t > 1
            # All non-depot pickup nodes at period t  are connected to the depot delivery node at period t+1
            A[2:end, 1, t, t+1, 1, 2] .= 1

            # All non-depot delivery nodes at period t  are connected to the depot pickup node at period t+1
            A[2:end, 1, t, t+1, 2, 1] .= 1

            # The depot delivery node at period t is connected to the depot delivery node at period t+1
            A[1, 1, t, t+1, 2, 2] = 1
        end


        if (t > 1) & (t < num_periods-1)
            # All non-depot pickup nodes at period t are connected to all non-depot pickup nodes at period t+1
            A[2:end, 2:end, t, t+1, 1, 1] .= 1

            # All non-depot delivery nodes at period t are connected to all non-depot nodes at period t+1
            A[2:end, 2:end, t, t+1, 2, :] .= 1
        end

        for i=1:num_nodes
            #  Delivery and pickup nodes at the same location and time are connected
            A[i, i, t, t, 2, 1] = 1
        end

    end
    return A
end


# Make cost array from pairwise node costs
function make_cost_array(
    num_nodes::Int64,
    num_periods::Int64,
    distance_matrix::Array{Int64,2}
)::Array{Int64, 6}
    C = zeros(Int, num_nodes, num_nodes, num_periods, num_periods, 2, 2)
    for i=1:num_nodes, j=1:num_nodes
        C[i, j, :, :, :, :] .= cost_matrix[i, j]
    end
    return C
end


# Solve network flow model
function solve_network_flow(
    node_balances::Array{Int64, 4},
    cost_matrix::Array{Int64, 2},
    num_locomotives::Int64,
    max_capacity::Int64;
    λ::Int64 = 10, # Penalty for longer tours
    log::Bool = true,
    solve_time::Int64 = 120
)::Tuple{Array{Int64,6}, Array{Int64,6}}

    if log
        println("Optimizing...")
        start_time = now()
    end

    # Process inputs
    num_nodes, num_periods, num_commodities, = size(node_balances)
    A = make_adjacency_array(num_nodes, num_periods)
    C = make_cost_array(num_nodes, num_periods, cost_matrix)
    total_supply = sum(node_balances[2:end, :, :, 1])
    total_demand = -sum(node_balances[2:end, :, :, 2])

    # Initialize model
    model = Model(solver=GurobiSolver(OutputFlag=0, TimeLimit=solve_time, GUROBI_ENV))

    # Define variables
    @variable(model, x[1:num_nodes, 1:num_nodes, 1:num_periods, 1:num_periods, 1:2, 1:2] >= 0)
    @variable(model, y[1:num_nodes, 1:num_nodes, 1:num_periods, 1:num_periods, 1:num_commodities, 1:2, 1:2] >= 0)
    @variable(model, z[1:num_nodes, 1:num_nodes, 1:num_periods, 1:num_periods, 1:2, 1:2], Bin)

    # Set network capacity constraints
    @constraint(model, z .<= A)
    @constraint(model, x[2:num_nodes, :, :, :, :, :] .<= max_capacity*z[2:num_nodes, :, :, :, :, :])
    @constraint(model, x[:, 2:num_nodes, :, :, :, :] .<= max_capacity*z[:, 2:num_nodes, :, :, :, :])
    @constraint(model, x[1, 1, :, :, :, 1] .<= total_demand*z[1, 1, :, :, :, 1])
    @constraint(model, x[1, 1, :, :, :, 2] .<= total_supply*z[1, 1, :, :, :, 2])

    #Set locomative conservation constraints
    @constraint(model, [t=1:num_periods-1],
        sum(z[:, :, t, t+1, :, :]) <= num_locomotives)
    @constraint(model, [i=2:num_nodes, t=1:num_periods-1, a=1:2],
        sum(z[:, i, t, t+1, :, :]) <= 1)
    @constraint(model, [i=2:num_nodes, t=1:num_periods, a=1:2],
        sum(z[:, i, :, t, :, a]) == sum(z[i, :, t, :, a, :]))

    # Set flow conservation constraints
    @constraint(model, [i=1:num_nodes, j=1:num_nodes, t=1:num_periods, u=1:num_periods, a=1:2, b=1:2],
        sum(y[i, j, t, u, :, a, b]) == x[i, j, t, u, a, b])
    @constraint(model, [i=1:num_nodes, t=1:num_periods,  k=1:num_commodities, a=1:2],
        sum(y[i, :, t, :, k, a, :]) - sum(y[:, i, :, t, k, :, a]) == node_balances[i, t, k, a])
    @constraint(model, [i=1:num_nodes, t=1:num_periods], sum(x[i, i, t, t, :, :]) == 0)

    # Set objective
    @objective(model, Min, sum(C.*z) + λ*sum(z[2:end, 2:end, :, :, :, :]))

    # Solve model
    solve(model, suppress_warnings=false)

    # Log performance if desired
    if log
        stop_time = now()
        objective_value = getobjectivevalue(model)
        objective_bound = getobjectivebound(model)
        println("Objective value: ", objective_value)
        println("Objective bound: ", objective_bound)
        println("MIP gap: ", round(
            (objective_value - objective_bound)/objective_bound * 100; digits=2)
        )
        println("Time elapsed: ", round(
            (stop_time - start_time).value/1e3; digits=2
        ), "s")
    end

    # Return raw solution
    return round.(Int, getvalue(x)), round.(Int, getvalue(z))
end


# Unpack raw solution from network flow solver into routes
function get_routes_from_raw_solution(
    x::Array{Int64, 6},
    z::Array{Int64, 6},
    node_balances::Array{Int64, 4},
    num_locomotives::Int64
)::Array{Array{Dict{Symbol,Int64},1},1}
    # Initialize data structures
    num_periods = size(node_balances, 2)
    route_nodes = zeros(Int, num_locomotives, num_periods)
    route_nodes[:, 1] .= 1
    route_nodes[:, end] .= 1
    route_status = ones(Int, num_locomotives, num_periods)
    route_deliveries = zeros(Int, num_locomotives, num_periods)
    route_pickups = zeros(Int, num_locomotives, num_periods)
    route_delivery_commodities = zeros(Int, num_locomotives, num_periods)
    route_pickup_commodities = zeros(Int, num_locomotives, num_periods)

    # For each period, infer the node, status and pickup delivery amounts for eahc route based on the active arcs
    for t=2:num_periods-1
        updated_routes = zeros(num_locomotives)
        for (from_node, to_node, from_activity, to_activity) in map(
            Tuple, findall(z[:, :, t-1, t, :, :] .== 1))
            # Find route corresponding to active arc
            if from_node == 1
                route_id = findall(
                    (route_nodes[:, t-1] .== from_node) .& (updated_routes .== 0)
                )[1]
            else
                route_id = findall(
                    (route_nodes[:, t-1] .== from_node) .& (route_status[:, t-1] .== from_activity)
                )[1]
            end

            # Update route node, status and delivery/pickup amounts
            route_nodes[route_id, t] = to_node
            if to_activity == 1
                if to_node != 1
                    route_pickups[route_id, t] = sum(node_balances[to_node, t, :, 1])
                    if route_pickups[route_id, t] != 0
                        route_pickup_commodities[route_id, t] = findall(
                            node_balances[to_node, t, :, 1] .> 0
                        )[1]
                    end
                end
                route_status[route_id, t] = 1
            else
                if to_node != 1
                    route_deliveries[route_id, t] = -sum(node_balances[to_node, t, :, 2])
                    if route_deliveries[route_id, t] != 0
                        route_delivery_commodities[route_id, t] = findall(
                            node_balances[to_node, t, :, 2] .< 0
                        )[1]
                    end
                end
                route_status[route_id, t] = 2
                if z[to_node, to_node, t, t, 2, 1] == 1
                    if to_node != 1
                        route_pickups[route_id, t] = sum(node_balances[to_node, t, :, 1])
                        if route_pickups[route_id, t] != 0
                            route_pickup_commodities[route_id, t] = findall(
                                node_balances[to_node, t, :, 1] .> 0
                            )[1]
                        end
                    end
                    route_status[route_id, t] = 1
                end
            end

            # Mark route as updated
            updated_routes[route_id] = 1
        end

        # Update other routes
        route_nodes[updated_routes .== 0, t] .= 1
    end

    # Return routes
    return [
        [
            Dict(
                :timestep => t,
                :node => route_nodes[route_id, t],
                :delivery_amount => route_deliveries[route_id, t],
                :pickup_amount => route_pickups[route_id, t],
                :delivery_commodity => route_delivery_commodities[route_id, t],
                :pickup_commodity => route_pickup_commodities[route_id, t]
            ) for t=1:num_periods
        ]
        for route_id in 1:num_locomotives
    ]
end


# Perform validation checks on routes (to be completed)
function validate_routes(
    routes::Array{Array{Dict{Symbol,Int64},1},1},
    node_balances::Array{Int64, 4}
)::Bool

    total_pickup_amount = sum([sum([stop[:pickup_amount] for stop in route]) for route in routes])
    total_delivery_amount = sum([sum([stop[:delivery_amount] for stop in route]) for route in routes])
    total_supply = sum(node_balances[2:end, :, :, 1])
    total_demand = -sum(node_balances[2:end, :, :, 2])

    if total_pickup_amount != total_supply
        println("Validation checks failed - pickup imbalance!")
        return false
    elseif total_delivery_amount != total_demand
        println("Validation checks failed - delivery imbalance!")
        return false
    else
        println("All validation checks passed.")
        return true
    end
end


# Optimize, unpack and validate routes
function optimize_routes(
    node_balances::Array{Int64, 4},
    cost_matrix::Array{Int64, 2},
    num_locomotives::Int64,
    max_capacity::Int64;
    λ::Int64 = 1,
    log::Bool = true,
    solve_time::Int64 = 120
)::Array{Array{Dict{Symbol,Int64},1},1}
    x_opt, z_opt = solve_network_flow(node_balances, cost_matrix, num_locomotives, max_capacity,
        λ=λ, log=log, solve_time=solve_time)
    routes = get_routes_from_raw_solution(x_opt, z_opt, node_balances, num_locomotives)
    @assert(validate_routes(routes, node_balances))
    return routes
end


# Process routes to produce data structures required for stage 2 optimization
function process_routes_for_stage_2(
    routes::Array{Array{Dict{Symbol,Int64},1},1}
)::Array{Tuple{Array{Int64, 1}, Array{Int64, 1}}, 1}
    num_periods = size(routes[1], 1)
    nodes = [[route[t][:node] for route in routes] for t=1:num_periods]
    return [(nodes[t], nodes[t+1]) for t=1:num_periods-1]
end
