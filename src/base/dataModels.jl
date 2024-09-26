using DataStructures
using NearestNeighbors
using StatsBase

mutable struct DataModel{T}
  data::Vector{T}
  func::Function
  DataModel{T}(func) where {T} = new{T}(Vector{T}(), func)
end

resetModels!(dm) = for symb in keys(dm)
    empty!(dm[symb].data)
  end

function emptyModels(; animation=false)
  dataModels = Dict{Symbol,DataModel}()
  if animation
    typ = Tuple{Float32,Vector{Int32},Vector{Int32},Vector{Int32}}
    dataModels[:animation] = DataModel{typ}(growAnimation)
  end
  return dataModels
end

function predefinedModels(dms=missing)::Dict{Symbol,DataModel}
  if ismissing(dms)
    dataModels = Dict{Symbol,DataModel}()
  else
    dataModels = dms
  end
  dataModels[:sectors] = DataModel{NTuple{4,Float64}}(sectorCoarsening)
  dataModels[:front] = DataModel{NTuple{5,Float64}}(frontMotion)

  return dataModels
end

"""
generate list of pairs -> (symbol, data)
"""
function unpack(dataModels)
  p = []
  foreach(x -> push!(p, Pair(x[1], x[2].data)), (dataModels))
  return p
end

"""
measures the mean and variance of the front height
"""
function frontMotion(active, graph, time)
  front = vcat(active[1].data, active[2].data, active[3].data)
  front_y = getfield.(graph[front], :y)
  m, v = mean_and_var(front_y)
  return (time, m, v, extrema(front_y)...)
end

"""
measures the mean and variance of the sector size at the front, and number of surviving roots
"""
function sectorCoarsening(active, graph, time)
  front = vcat(active[1].data, active[2].data, active[3].data)
  uniqueElems = countmap(getfield.(graph[front], :ID_2))
  nroots = length(keys(uniqueElems))
  m, v = mean_and_var(values(uniqueElems))
  return (time, m, v, nroots)
end

"""
take snapshots of the front to later make an animation
"""
function growAnimation(active, graph, time)
  lineages = Int32[]
  if time > 0.0
    front = vcat(active[1].data, active[2].data, active[3].data, active[4].data)
    phylogeny, _ = genealogy(graph, front)
    lineages = collect(keys(phylogeny))
  end
  ID_1 = getfield.(graph, :ID_2)
  ID_2 = getfield.(graph, :ID_3)
  return (time, ID_1, ID_2, lineages)
end

# using P.B.C get the arc distance between two points
""" 
enforce periodic boundary conditions along X
"""
function periodicBoundaries(halfwidth::Int, width::Int)
  f = let h = halfwidth, w = width
    x -> (x > h) ? w - x : x
  end
  return f
end

"""
construct ancestry trees as a dictionary (julia implemented as hash table)
"""
function genealogy(graph, front)
  phylogeny = Dict{Int64,Int64}()

  # associate originID to list of branch splits
  branching = DefaultOrderedDict{Int64,Vector{Int}}([])

  @inbounds for nodeIdx in front
    ID = graph[nodeIdx].ID_2
    current = nodeIdx
    while true
      if haskey(phylogeny, current)
        push!(branching[ID], current)
        break
      end

      phylogeny[current] = graph[current].ancestor
      current = graph[current].ancestor
      current == 0 && break
    end
  end

  foreach(unique!, values(branching))
  return phylogeny, branching
end

"""
provided with a list of position tuples (input), find all indices of targets for which an element is within a distance radius of an input element.
"""
function overlap(radius, xy, srcIdxs, kdtree)
  length(kdtree.data) == 0 && return zeros(Int64, 0)

  # created containers for tracking which  targets get visited
  reached = zeros(Bool, length(kdtree.data))

  # better to query with all sample points
  # task treat edge cases, periodic boundary conditions
  # -> (idx in kdtree.data, element in kdtree.data w/ min dist )
  idxs, dist = nn(kdtree, SVector.(xy[srcIdxs]))

  # update visited targets by input
  @inbounds foreach(i -> reached[i] = true, idxs[dist .<= radius])
  return reached
end
