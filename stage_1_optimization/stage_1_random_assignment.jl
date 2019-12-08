using Random

function get_feasible_locomotives_for_job(job, locomotives, max_capacity, max_consecutive_deliveries=3)
    feasible_locomotives = []
    for locomotive in locomotives
        if job[:delivery_amount] > 0
            if (locomotive[:status] == 0) & (job[:delivery_amount] + locomotive[:cars] < max_capacity) & (
                locomotive[:delivery_count] < max_consecutive_deliveries)
                push!(feasible_locomotives, locomotive[:id])
            end
        else
            if (locomotive[:status] == 0) | ((locomotive[:status] == 1) & (
                        job[:pickup_amount] + locomotive[:cars] < max_capacity))
                push!(feasible_locomotives, locomotive[:id])
            end
        end
    end
    return feasible_locomotives
end


function assign_locomotive_to_job(job, locomotives, max_capacity)
    feasible_locomotives = get_feasible_locomotives_for_job(job, locomotives, max_capacity)
    return rand(feasible_locomotives)
end

function update_locomotive_for_job(job, locomotive)
    if job[:pickup_amount] > 0
        new_status = 1
        new_delivery_count = 0
    else
        new_status = 0
        new_delivery_count = locomotive[:delivery_count] + 1
    end

    if (new_status ==  1) & (locomotive[:status] == 1)
        new_cars = locomotive[:cars] + job[:pickup_amount]
    elseif (new_status ==  1) & (locomotive[:status] == 0)
        new_cars = job[:pickup_amount]
    else
        new_cars = locomotive[:cars] + job[:delivery_amount]
    end

    return Dict(
        :id => locomotive[:id],
        :status => new_status,
        :cars => new_cars,
        :route => vcat(locomotive[:route], job[:node]),
        :delivery_count => new_delivery_count
    )
end


function update_unassigned_locomotive(locomotive, t, num_periods)
    if t < num_periods - 1
        new_node = rand([1, locomotive[:route][end]])
    else
        new_node = 1
    end
    if new_node == 1
        new_status = 0
        new_cars = 0
        new_delivery_count = 0
    else
        new_status = locomotive[:status]
        new_cars = locomotive[:cars]
        new_delivery_count = locomotive[:delivery_count] + 1
    end
    return Dict(
        :id => locomotive[:id],
        :status => new_status,
        :cars => new_cars,
        :route => vcat(locomotive[:route], new_node),
        :delivery_count => new_delivery_count
    )
end


function get_jobs_from_routes(routes)
    jobs = vcat(
        [
            [stop for stop in route if max(stop[:pickup_amount], stop[:delivery_amount]) > 0]
            for route in routes
        ]...
    )
    jobs = sort(collect(jobs), by=job->job[:timestep])
    return jobs
end


function initialize_locomotives(num_locomotives)
    return [
        Dict(
            :id => id,
            :status => 0,
            :cars => 0,
            :route => [1],
            :delivery_count => 0
        )
        for id=1:num_locomotives
    ]
end


function generate_random_assigments(jobs, num_periods, max_capacity)
    locomotives = initialize_locomotives(num_locomotives)
    for t=1:num_periods
        assigned_locomotive_indicator = zeros(num_locomotives)
        for job in [job for job in jobs if (job[:timestep] == t)]
            selected_locomitive_index = assign_locomotive_to_job(job,
                            locomotives[findall(assigned_locomotive_indicator .== 0)], max_capacity)
            locomotives[selected_locomitive_index] = update_locomotive_for_job(job,
                            locomotives[selected_locomitive_index])
            assigned_locomotive_indicator[selected_locomitive_index] = 1
        end
        for l in findall(assigned_locomotive_indicator .== 0)
            locomotives[l] = update_unassigned_locomotive(locomotives[l], t, num_periods)
        end
    end
    return locomotives
end


function get_cost_for_optimal_routes(routes, cost_matrix, num_periods)
    cost = 0
    for route in routes
        for t=1:num_periods-1
            from = route[t][:node]
            to =  route[t+1][:node]
            cost += cost_matrix[from, to]
        end
    end
    return cost
end


function get_total_cost_of_random_assignments(locomotives, cost_matrix, num_periods)
    cost = 0
    for locomotive in locomotives
        for t=2:num_periods
            from = locomotive[:route][t]
            to = locomotive[:route][t+1]
            cost += cost_matrix[from, to]
        end
    end
    return cost
end


function generate_random_assignment_costs(jobs, num_periods, max_capacity, num_trials=10000)
    random_assignment_costs = []
    for trial=1:num_trials
        try
            locomotives = generate_random_assigments(jobs, num_periods, max_capacity)
            push!(
                random_assignment_costs,
                get_total_cost_of_random_assignments(locomotives, cost_matrix, num_periods)
            )
        catch
        end
    end
    return random_assignment_costs
end


function estimate_pct_cost_saving(min_cost, random_assignment_costs)
    mean_cost = sum(random_assignment_costs)/length(random_assignment_costs)
    return (mean_cost - min_cost)/mean_cost * 100
end

function plot_cost_distribution(min_cost, random_assignment_costs)
    plt.axvline(min_cost, color="tab:orange", label="Min cost")
    plt.hist(random_assignment_costs, density=true)
    plt.xlabel("Cost")
    plt.ylabel("Normalized frequency");
end
