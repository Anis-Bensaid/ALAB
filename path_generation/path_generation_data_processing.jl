using CSV, DataFrames, Suppressor, LinearAlgebra
#include("path_generation.jl")

# Compute cost matrix based on shortest paths for specified terminals
function compute_cost_matrix_using_shortest_path(
    NodeCoords::Array{Any,2},
    Edges::Array{Int64,2},
    terminalNames::Array{String, 1}
)::Array{Int64,2}
    distMat = EuclDist(NodeCoords[:,2:3], NodeCoords[:,2:3]).*AdjMat(Edges, size(NodeCoords,1))
    terminalNodes = findall(NodeCoords[:,4] .== "TRUE")
    numTerminals = size(terminalNodes, 1)
    costMat = zeros(Int, numTerminals, numTerminals)
    allPaths = 0
    for i=1:numTerminals
        allPaths  = @suppress K_ShortestPaths(
            1, NodeCoords, Edges, distMat, repeat([i], numTerminals), terminalNodes, false);
        costMat[i,:] = Int.(round.([allPaths[j][1][end].f for j=1:numTerminals]))
    end
    costDf = DataFrame(costMat);
    costDf[!,:from] =  NodeCoords[terminalNodes,1]
    numIncludedTerminals = size(terminalNames, 1)
    orderedSubsetCostMat = zeros(Int, numIncludedTerminals, numIncludedTerminals)
    for i=1:numIncludedTerminals
        for j=1:numIncludedTerminals
            orderedSubsetCostMat[i,j] = costDf[
                costDf[:,:from] .== terminalNames[i],
                findall(costDf[:, :from] .== terminalNames[j])[1]
            ][1]
        end
    end
    return Int.(round.((orderedSubsetCostMat .+ orderedSubsetCostMat')/2))
end

# Compute cost matrix based on euclidian between specified terminals
function compute_cost_matrix_using_euclidian_distance(
    NodeCoords::Array{Any,2},
    terminalNames::Array{String, 1}
)::Array{Int64,2}
    numIncludedTerminals = size(terminalNames, 1)
    costMat = zeros(Int, numIncludedTerminals, numIncludedTerminals)
    for i=1:numIncludedTerminals
        for j=1:numIncludedTerminals
            from = NodeCoords[NodeCoords[:,1].== terminalNames[i], 2:3]
            to = NodeCoords[NodeCoords[:,1].== terminalNames[j], 2:3]
            costMat[i,j] = Int(round(norm(to .- from)))
        end
    end
    return costMat
end

# Get terminal names and commodities from schedule
function get_terminal_and_commodity_names_from_schedule(
    deliveries_file_path::String
)::Tuple{Array{String, 1}, Array{String, 1}}
    deliveries = CSV.read(deliveries_file_path)
    terminal_names = unique(deliveries.Location)
    commodity_names = unique(deliveries.Commodity)
    return terminal_names, commodity_names
end


# Get cost matrix for corresponding to terminals in schedule
function get_cost_matrix(
    nodes_file_path::String,
    edges_file_path::String,
    deliveries_file_path::String,
    method::String = "euclidian"
)::Array{Int64,2}
    NodeCoords = convert(Matrix, CSV.read(nodes_file_path, header = true))
    Edges = convert(Matrix, CSV.read(edges_file_path, header = true))[:,2:3]
    terminalNames, = get_terminal_and_commodity_names_from_schedule(deliveries_file_path)
    if method == "euclidian"
        return compute_cost_matrix_using_euclidian_distance(NodeCoords, terminalNames)
    elseif method == "shortest path"
        return compute_cost_matrix_using_shortest_path(NodeCoords, Edges, terminalNames)
    else
        @assert(false, "Method not implemented.")
    end
end
