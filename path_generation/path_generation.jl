# RUMO LOGISTICA ROUTING PROBLEM (K-SHORTEST PATHS)

# Import required packages -----------------------------------------------------
using JuMP, Gurobi, LinearAlgebra, DataStructures;
using Plots; theme(:wong2);
#clibrary(:colorbrewer);

## ALGORITHM IMPLEMENTATION ####################################################

# Define the properties of the nodes -------------------------------------------
mutable struct Node
    name::String
    g::Float64
    h::Float64
    f::Float64
    Parent::Union{Node, Missing}
    Position::Tuple{Float64, Float64}
end
name(N::Node) = N.name
g(N::Node) = N.g
h(N::Node) = N.h
f(N::Node) = N.f
Parent(N::Node) = N.Parent
Position(N::Node) = N.Position

equal(P1::Node, P2::Node) = Position(P1) == Position(P2)
updateCost(N::Node, val, string) = if isequal(string, "g") N.g = val
                                   elseif isequal(string, "h") N.h = val
                                   elseif isequal(string, "f") N.f = val end
updateParent(P1::Node, P2::Node) = P1.Parent = P2

# Define the properties of the edges -------------------------------------------
mutable struct Edge
    name::String
    source::Node
    sink::Node
    dist::Float64
end
name(E::Edge) = E.name
source(E::Edge) = E.source
sink(E::Edge) = E.sink
dist(E::Edge) = E.dist

# Define the properties of the graph -------------------------------------------
mutable struct Graph
    nodes::Array{Union{Any, Node}, 1}
    edges::Dict{}
end
nodes(G::Graph) = G.nodes
edges(G::Graph) = G.edges
addNode(G::Graph, node::Node) = if true push!(G.nodes, node)
                                if true G.edges[node.name] = [] end end
addEdge(G::Graph, edge::Edge) = push!(G.edges[edge.source.name], edge)
removeEdge(G::Graph, P1::Node, P2::Node) = filter!(x->x.sink.name!= P2.name, G.edges[P1.name])
removeNode(G::Graph, P1::Node) = delete!(G.edges, P1.name)
childrenOf(G::Graph, node::Node) = G.edges[node.name]

# ASTAR algorithm -----------------------------------------------------------
function AStar(network, startNode, endNode)
    # Returns a list of Nodes as path from start to final node

    # Initialize open and closed lists
    openList = PriorityQueue()
    closedList = []

    if isempty(childrenOf(network, startNode)) || isempty(childrenOf(network, endNode))
        println("The starting or terminal stations are not in the newtork")
        return missing, Inf
    end

    # Add start node
    enqueue!(openList, startNode, startNode.f)

    # Loop until reach sink node
    while !isempty(openList)

        currentNode = dequeue!(openList)
        push!(closedList, currentNode)

        # If at destination, return path
        if currentNode.name == endNode.name
            path = []
            current = currentNode
            while current !== missing
                push!(path, current)
                current = current.Parent
            end
            return(reverse!(path), path[end].f)
        end

        # Look at neighbors of current node
        for neighbor in childrenOf(network, currentNode)

            for closed_child in closedList
                if equal(neighbor.sink, closed_child)
                    continue
                end
            end

            temp_score = currentNode.g + neighbor.dist

            if temp_score < neighbor.sink.g
                for open in openList
                    if equal(neighbor.sink, open[1])
                        continue
                    end
                end
                try
                    neighbor.sink.g = temp_score
                    val = sqrt((neighbor.sink.Position[1] - endNode.Position[1])^2 +
                          (neighbor.sink.Position[2] - endNode.Position[2])^2)
                    updateCost(neighbor.sink, val, "h")
                    updateCost(neighbor.sink, neighbor.sink.g + neighbor.sink.h, "f")
                    updateParent(neighbor.sink, currentNode)
                    enqueue!(openList, neighbor.sink => neighbor.sink.f)
                catch
                    @warn "Node not added"
                end
            end
        end
    end
    return (missing, Inf)
end

## DEFINE HELPER FUNCTIONS #####################################################
function AdjMat(edges, n)
    adj = zeros(n,n)
    for k=1:size(edges, 1)
        adj[edges[k,1], edges[k,2]] = 1
    end
    return Int.(adj)
end

function EuclDist(A, B)
 # Determines the Eucledian Distance between two matrices
    (m, n) = size(A)
    (p, q) = size(B)
    dist = zeros(m, p)
    for i in 1:size(A,1)
        for j in 1:size(B,1)
          dist[i, j] = sqrt(sum((A[i, :].-B[j,:]).^2))
        end
    end
    return(dist)
end

function BuildGraph(coord, branches, distMat)
    graph = Graph([], Dict())
    for i = 1:size(coord, 1)
        node = Node(coord[i,1], Inf,0, Inf, missing, (coord[i,2], coord[i,3]))
        addNode(graph, node)
    end

    for j = 1:size(branches, 1)
        ic_1 = branches[j,1]
        ic_2 = branches[j,2]
        edge = Edge(string(j), graph.nodes[ic_1], graph.nodes[ic_2], distMat[ic_1, ic_2])
        rev_edge = Edge(string(j), graph.nodes[ic_2], graph.nodes[ic_1], distMat[ic_1, ic_2])
        addEdge(graph, edge)
        addEdge(graph, rev_edge)
    end
    return graph
end

function PlotInitialGraph(NodeCoords, Edges)
    P = scatter(NodeCoords[:,2], NodeCoords[:,3])
    for b=1:size(Edges,1)
        #x = [NodeCoords[findfirst(x->x==i, NodeCoords[:,1]), 2] for i in Edges[b, :]]
        #y = [NodeCoords[findfirst(x->x==i, NodeCoords[:,1]), 3] for i in Edges[b, :]]
        x = [NodeCoords[i,2] for i in Edges[b, :]]
        y = [NodeCoords[i,3] for i in Edges[b, :]]
        plot!(P, x, y, color = "black", linewidth = 0.6, legend = false)
    end
    return P
end

function plotPath(P, nodes::Vector{Any})
    x = [node.Position[1] for node in nodes]
    y = [node.Position[2] for node in nodes]
    plot!(P, x, y, alpha = 0.4, linewidth = 3)
    return P
end

function K_ShortestPaths(K, NodeCoords, Edges, distMat, sources, sinks, toPrint)
    K::Int64
    NodeCoords::Array{Any, 2}
    Edges::Array{Int64, 2}
    distMat::Array{Float64, 2}
    sources::Array{Int64}
    sinks::Array{Int64}
    toPrint::Bool

    # Check that matrix dimensions match
    if size(distMat,1) != size(distMat,2) || size(NodeCoords, 1) != size(distMat,2)
        error("Check dimensions of distance matrix")
    end

    graph = BuildGraph(NodeCoords, Edges, distMat)
    allPaths = []

    for i = 1:length(sources)
        KBest = []
        Potential = PriorityQueue()

        # Make a copy of graph and initialize nodes
        graphCopy = deepcopy(graph)
        Nodes = nodes(graphCopy)
        startNode = Nodes[sources[i]]; startNode.g = 0; startNode.f = 0;
        endNode = Nodes[sinks[i]]

        # Run the A* algorithm
        (shortest, value) = AStar(graphCopy, startNode, endNode)
        push!(KBest, shortest)

        # Begin the Yen's KSP algorithm
        for k = 2:K
            for j = 1:size(KBest[k-1],1)-1
                tempGraph = deepcopy(graph)
                spurNode = KBest[k-1][j]
                rootPath = KBest[k-1][1:j]

                for path in KBest
                        if length(path)>length(rootPath) && rootPath == path[1:j]
                            removeEdge(tempGraph, path[j], path[j+1])
                            removeEdge(tempGraph, path[j+1], path[j])
                        end
                end

                for idx = 1:length(rootPath)-1
                    removeEdge(tempGraph, rootPath[idx], rootPath[idx+1])
                    removeEdge(tempGraph, rootPath[idx+1], rootPath[idx])
                end


                newPath, val = AStar(tempGraph, spurNode, endNode)
                if newPath !== missing
                    totPath = vcat(rootPath[1:end-1], newPath)
                    enqueue!(Potential, totPath, val + rootPath[end].f)
                end
            end
            if isempty(Potential)
                break
            end
            push!(KBest, dequeue!(Potential))
        end
        push!(allPaths, KBest)
    end

    if toPrint
        for routes in allPaths
            println("------------------------------")
            println("Begin at position: ", routes[1][1].Position, ", End at position: ", routes[1][end].Position)
            println("\nThe ordered shortest paths are\n ")
            for path in routes
                for N in path
                    print(N.name, " => ")
                end
                print(" Distance travelled: ", path[end].f, "\n")
            end
        end
    end
    return allPaths
end
