# doc: IDs(environment condition | ancestral lineage | mutant type | mutant number)
mutable struct Node
  ancestor::UInt32
  filled::Bool
  ID_1::Int32
  ID_2::Int32
  ID_3::Int32
  ID_4::UInt32
  time::Float32
  nbors::UInt32
  Node(n) = new(0, false, 0, 0, 0, 0, 0.0, n)
end

@inline function updateNeighbors!(graph, cntns, nodeIndx, active, container)
  for nbor in view(cntns, :, nodeIndx)
    nbor == 0 && continue

    # subtract one from all neighbors
    graph[nbor].nbors -= 1

    # if this site, or neighbor, is surrounded, then remove it
    if graph[nbor].nbors == 0 && graph[nbor].filled
      HashMaps.remove!(active[graph[nbor].ID_1], nbor)
      container[graph[nbor].ID_1] -= 1
    end
  end
end

resetGraph!(graph, cntns) = @inbounds @simd for i in eachindex(graph)
  graph[i].filled = false
  graph[i].ancestor = 0
  graph[i].ID_1 = 0
  graph[i].ID_2 = 0
  graph[i].ID_3 = 0
  graph[i].ID_4 = 0
  graph[i].time = 0.0
  graph[i].nbors = sum(view(cntns, :, i) .> 0)
end

gNodes(typ::Symbol=:hex) =
  if typ == :hex
    # https://www.redblobgames.com/grids/hexagons/#neighbors
    bottom = [(1, 0), (-1, 0), (-1, -1), (0, -1), (-1, 1), (0, 1)]
    top = [(1, 0), (-1, 0), (0, 1), (1, 1), (0, -1), (1, -1)]
    return Dict{Int,Vector{NTuple{2,Int}}}(1 => bottom, 0 => top)
  end

function hex_grid_connections()
  # https://www.redblobgames.com/grids/hexagons/#neighbors
  bottom = [(1, 0), (-1, 0), (-1, -1), (0, -1), (-1, 1), (0, 1)]
  top = [(1, 0), (-1, 0), (0, 1), (1, 1), (0, -1), (1, -1)]
  return Dict{Int,Vector{NTuple{2,Int}}}(1 => bottom, 0 => top)
end

function square_grid_connections()
  bottom = [(-1, 0), (1, 0), (0, -1), (0, 1)]
  top = [(-1, 0), (1, 0), (0, -1), (0, 1)]
  return Dict{Int,Vector{NTuple{2,Int}}}(1 => bottom, 0 => top)
end

function constructGraph(dimensions, nodes, T::Type, cnstr, periodic=true)
  lx, ly = dimensions
  graph = Vector{T}(undef, lx * ly)
  connections = zeros(UInt32, 6, lx * ly)

  for row in 1:ly, col in 1:lx
    nodeIndx = lx * (row - 1) + col
    nbors::UInt32 = 0

    for (i, (dx, dy)) in enumerate(nodes[row % 2])
      ny = row + dy
      0 < ny <= ly || continue

      if periodic
        nx = mod(col + dx, 1:lx)
      else
        nx = col + dx
        0 < nx <= lx || continue
      end

      idx = lx * (ny - 1) + nx
      connections[i, nodeIndx] = idx
      nbors += 1
    end
    graph[nodeIndx] = cnstr(nbors, col, row)
  end
  return graph, connections
end

hexGraph(dims, periodic=true) =
  let
    cvr = (x, y) -> (sqrt(3) * (x - 0.5 * (y % 2)), 1 + 1.5 * (y - 1))
    constructGraph(dims, hex_grid_connections(), Node, (n, x, y) -> Node(n), periodic)
  end

squareGraph(dims, periodic=true) =
  let
    constructGraph(dims, square_grid_connections(), Node, (n, x, y) -> Node(n), periodic)
  end
