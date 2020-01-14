using CSV, DataFrames

function get_node_balances_from_schedule(
    deliveries_file_path::String,
    pickups_file_path::String
)::Array{Int64, 4}
    # Load CSV files
    deliveries = CSV.read(deliveries_file_path)
    pickups = CSV.read(pickups_file_path)

    # Get problem dimensions
    num_nodes = length(unique(deliveries.Location))
    num_periods = size(deliveries, 2) - 2
    num_commodities = length(unique(deliveries.Commodity))

    # Set node balances
    node_balances = zeros(Int, num_nodes, num_periods, num_commodities, 2)
    for (i, location) in enumerate(groupby(deliveries, :Location))
        if i == 1
            node_balances[i, :, :, 1] = convert(Matrix, location)[:,3:end]'
        else
            node_balances[i, :, :, 2] = -convert(Matrix, location)[:,3:end]'
        end
    end
    for (i, location) in enumerate(groupby(pickups, :Location))
        if i == 1
            node_balances[i, :, :, 2] = -convert(Matrix, location)[:,3:end]'
        else
            node_balances[i, :, :, 1] = convert(Matrix, location)[:,3:end]'
        end
    end
    @assert(sum(node_balances) == 0)
    return node_balances
end
