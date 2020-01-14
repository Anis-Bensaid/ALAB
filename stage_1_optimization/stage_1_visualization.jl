using PyPlot

# Print sequential steps of route (a bit hacky)
function print_route(
    route::Array{Dict{Symbol,Int64},1},
    num_periods::Int64,
    terminal_names::Array{String, 1},
    commodity_names::Array{String, 1};
    start_hour::Int64 = 8
)
    spaces = 10
    println(repeat('=', 100))
    time_strings = [string(0, hour,":00")[end-4:end] for hour=start_hour:start_hour+num_periods]
    for stop in route
        if stop[:node] == 1
            println(
                "Time: ", time_strings[stop[:timestep]],
                "\t Terminal: ", string(terminal_names[1], repeat(" ", spaces))[1:spaces],
                "\t Activity: load/unload"
            )
        elseif (stop[:delivery_amount] == 0) & (stop[:pickup_amount] == 0)
            println(
                "Time: ", time_strings[stop[:timestep]],
                "\t Terminal: ", string(terminal_names[stop[:node]], repeat(" ", spaces))[1:spaces],
                "\t Activity: wait"
            )
        else
            if stop[:delivery_amount] > 0
                println(
                    "Time: ", time_strings[stop[:timestep]],
                    "\t Terminal: ", string(terminal_names[stop[:node]], repeat(" ", spaces))[1:spaces],
                    "\t Activity: deliver ", stop[:delivery_amount], " cars of ",
                    commodity_names[stop[:delivery_commodity]],
                )
            end
            if stop[:pickup_amount] > 0
                println(
                    "Time: ", time_strings[stop[:timestep]],
                    "\t Terminal: ", string(terminal_names[stop[:node]], repeat(" ", spaces))[1:spaces],
                    "\t Activity: pick-up ", stop[:pickup_amount], " cars of ",
                    commodity_names[stop[:pickup_commodity]],
                )
            end
        end
    end
end


# Plot route assignment solution
function plot_routes(
    routes::Array{Array{Dict{Symbol,Int64},1},1},
    node_balances::Array{Int64, 4},
    terminal_names::Array{String, 1};
    start_hour::Int64 = 8,
    node_size::Int64 = 200,
    fig_size::Tuple{Int64, Int64} = (12, 8)
)::Figure

    # Get plotting inputs
    num_nodes, num_periods, = size(node_balances)
    num_locomotives = size(routes, 1)
    locations_index = vec(1:num_nodes) .- 1
    timesteps =start_hour:(start_hour + num_periods - 1)
    activity_info = [
        (1, "Pick-up", "whitesmoke"),
        (2, "Delivery", "lightgrey"),
        (3, "Delivery & pick-up", "dimgrey")
    ]
    activities = zeros(Int, num_nodes, num_periods)
    for i=1:num_nodes
        for t=1:num_periods
            delivery_indicator = sum(node_balances[i,t,:,2]) < 0
            pickup_indicator = sum(node_balances[i,t,:,1]) > 0
            if delivery_indicator & pickup_indicator
                activities[i, t] = 3
            elseif delivery_indicator
                activities[i, t] = 2
            elseif pickup_indicator
                activities[i, t] = 1
            end
        end
    end
    fig = plt.figure(figsize=fig_size)

    #Route leg line plots
    line_plots = []
    for l=1:num_locomotives
        line_plot, = plt.plot(
            timesteps,
            [locations_index[i] for i in [stop[:node] for stop in routes[l]]])
        push!(line_plots, line_plot)
    end

    # Activity scatter plots
    scatter_plots = []
    for t=1:num_periods
        scatter_plots = []
        for (activity_index, activity_name, activity_color) in activity_info
            selected_locations = locations_index[2:end][activities[2:end,t] .== activity_index]
            push!(
                scatter_plots,
                plt.scatter(
                    timesteps[t]*ones(size(selected_locations, 1)),
                    selected_locations,
                    color=activity_color,
                    s=node_size,
                    edgecolor="black",
                    zorder=10
                )
            )
        end
    end
    push!(
        scatter_plots,
        plt.scatter(timesteps, zeros(num_periods), color="black", s=node_size, zorder=10)
    )

    # Legends
    legend1 = plt.legend(
        line_plots,
        [string("Locomotive ", l) for l=1:num_locomotives],
        bbox_to_anchor=(1.175,1)
    )
    plt.legend(
        scatter_plots,
        vcat([activity_name for (activity_index, activity_name, activity_color) in activity_info], "KM5"),
        bbox_to_anchor=(1,0.175),
    )
    plt.gca().add_artist(legend1)

    # Axis labels
    plt.xlabel("Time", fontsize=14)
    plt.ylabel("Terminal", fontsize=14)
    plt.xticks(timesteps, [string(time,":00") for time in timesteps])
    plt.yticks(locations_index, terminal_names)
    return fig
end
