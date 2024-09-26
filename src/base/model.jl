module Model
using StatsBase
include("containers.jl")
import ArgParse: ArgParseSettings, parse_args, add_arg_table!, parse_item, add_arg_group!
import StaticArrays: SVector
import .HashMaps: HashVec, add!, remove!
export simulate!, resetGraph!, parseArgs, setPath, hexGraph

parse_item(::Type{NTuple{3,T}}, x::AbstractString) where {T} = Tuple(convert.(T, parse.(Float64, split(x, ','))))

function parseArgs()
  addArgs!(sts, name; kwargs...) = add_arg_table!(sts, "--$(name)", Dict(kwargs))

  sts = ArgParseSettings()
  addArgs!(sts, "landscape"; required=false, arg_type=String)
  addArgs!(sts, "env_type"; arg_type=String, default="uniform")
  addArgs!(sts, "separation"; arg_type=Int)
  addArgs!(sts, "gap"; arg_type=Int, default=0)

  addArgs!(sts, "numberTrials"; arg_type=Int, default=1)
  addArgs!(sts, "numberSamples"; arg_type=Int, default=50)

  addArgs!(sts, "width"; required=true, arg_type=Int64)
  addArgs!(sts, "height"; required=true, arg_type=Int64)
  addArgs!(sts, "selection"; required=true, arg_type=Float64)
  addArgs!(sts, "mutation"; required=true, arg_type=Float64)
  addArgs!(sts, "intensity"; required=true, arg_type=Float64)
  addArgs!(sts, "density"; required=false, arg_type=Float64, default=0.09)
  addArgs!(sts, "radius"; required=false, arg_type=Int64, default=10)

  addArgs!(sts, "outputPath"; required=true, arg_type=String)
  addArgs!(sts, "animate"; action=:store_true)
  addArgs!(sts, "rewrite"; action=:store_true)
  addArgs!(sts, "rngSeed"; arg_type=Int, default=1)
  addArgs!(sts, "printInfo"; action=:store_true)
  addArgs!(sts, "detailed_analytics"; action=:store_true)
  addArgs!(sts, "heatmap"; action=:store_true)
  addArgs!(sts, "initial_type"; arg_type=String, default="alt")

  # exclusive option groups
  add_arg_group!(sts, "exlusive"; exclusive=true)
  addArgs!(sts, "standing_variation"; action=:store_true)
  addArgs!(sts, "singleMutant"; action=:store_true)
  parsedArgs = namedtuple(parse_args(sts))

  if parsedArgs.printInfo
    println("Simulation Info:")
    foreach(k -> println("  > $k  =>  $(parsedArgs[k])"), keys(parsedArgs))
  end
  return parsedArgs
end

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

gNodes(typ::Symbol=:hex) =
  if typ == :hex
    # https://www.redblobgames.com/grids/hexagons/#neighbors
    bottom = [(1, 0), (-1, 0), (-1, -1), (0, -1), (-1, 1), (0, 1)]
    top = [(1, 0), (-1, 0), (0, 1), (1, 1), (0, -1), (1, -1)]
    return Dict{Int,Vector{NTuple{2,Int}}}(1 => bottom, 0 => top)
  end

graphSources(cli) = collect(1:(cli.width)) .+ cli.width * (cli.height - 1)

function populate!(graph, cntns, dims, env, opts; row=1, num=3, standingVar=false)
  active = Vector{HashVec{UInt32}}(undef, num)
  foreach(i -> active[i] = HashVec{UInt32}(), collect(1:1:num))

  if standingVar
    if contains(opts.initial_type, "split")
      strainID = ones(Int, dims[1])
      strainID[cld(dims[1], 2):end] .= 2
    elseif contains(opts.initial_type, "alt")
      strainID = mod.(collect(1:dims[1]), 2) .+ 1
    else
      strainID = rand([1, 2], dims[1])
    end
  else
    strainID = ones(Int, dims[1])
  end

  for col in 1:dims[1]
    # shift group affliation to match environment
    if strainID[col] == 1
      nodeID = env[col] == 1 ? 1 : 2
    else
      nodeID = env[col] == 1 ? 3 : 4
    end

    nodeIndx = dims[1] * (row - 1) + col

    graph[nodeIndx].filled = true
    graph[nodeIndx].ID_1 = nodeID
    graph[nodeIndx].ID_2 = col
    graph[nodeIndx].ID_3 = strainID[col]

    # don't add site to front if it already has no empty nbors
    if graph[nodeIndx].nbors > 0
      add!(active[nodeID], nodeIndx)
    end

    for nbor in view(cntns, :, nodeIndx)
      nbor == 0 && continue

      # subtract one from all neighbors
      graph[nbor].nbors -= 1

      # if this site, or neighbor, is surrounded, then remove it
      if graph[nbor].nbors == 0 && graph[nbor].filled
        remove!(active[graph[nbor].ID_1], nbor)
      end
    end
  end
  return active
end

function buildMap(dimensions, nodes, T::Type, cnstr, periodic=true)
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
    buildMap(dims, gNodes(), Node, (n, x, y) -> Node(n), periodic)
  end

@inline function nextIndex(graph, cntns, idx)
  nbors = UInt32[]
  @inbounds for nbor in view(cntns, :, idx)
    (nbor == 0 || graph[nbor].filled) && continue
    push!(nbors, nbor)
  end
  if length(nbors) == 1
    return first(nbors)
  elseif length(nbors) == 0
    error("found node with zero neighbors $idx")
  end
  return rand(nbors)
end

function ssa(active, rates)::Tuple{Int,Float64}
  a = length(active[1].data) * rates[1]
  b = length(active[2].data) * rates[2]
  c = length(active[3].data) * rates[3]
  r = a + b + c + length(active[4].data) * rates[4]
  η = rand() * r

  if η < a
    return 1, r
  elseif η < a + b
    return 2, r
  elseif η < a + b + c
    return 3, r
  else
    return 4, r
  end
end

function simulate!(graph, cntns, cli, dataModels, env)
  sel = cli.selection
  intensity = cli.intensity
  mut = cli.mutation

  rates = SVector{4,Float64}([1, 1 + intensity, 1 - sel, (1 - sel) * (1 + intensity)])

  width = cli.width
  height = cli.height

  time = 0.0

  active = populate!(graph, cntns, (width, height), env, cli; num=4, standingVar=cli.standing_variation)

  nrecord::Int64 = cld(width * height, cli.numberSamples)
  itrCntr::Int64 = nrecord
  continueRecording::Bool = true

  # .track number of mutations that occur
  mutationCount = 0

  # .save snapshot of the front as it touches the top
  front = Int32[]
  stopCondition = width * (height - 1) + 1

  tExtinction = -1

  willMutate = true
  singleMutant = cli.singleMutant || cli.env_type == "circle"

  @inbounds while true
    length(active[1].data) + length(active[2].data) + length(active[3].data) + length(active[4].data) == 0 && break

    if (length(active[1].data) + length(active[2].data)) == 0 && tExtinction < 0
      tExtinction = time
    end

    if continueRecording && itrCntr == nrecord
      itrCntr = 0
      for (_, measure) in dataModels
        push!(measure.data, measure.func(active, graph, time))
      end
    end
    itrCntr += 1

    # .find which cell will grow next
    groupID, R = ssa(active, rates)

    # .parent and child index
    parentIdx = rand(active[groupID].data)
    childIdx = nextIndex(graph, cntns, parentIdx)

    # .update time
    time += -log(rand()) / R

    # .mutate childID based on mutation rate, and increment mutations counts
    mutID = graph[parentIdx].ID_3
    if !singleMutant && mutID == 1 && rand() < mut
      mutID = 2
      mutationCount += 1
    end

    if singleMutant && willMutate && graph[parentIdx].ID_3 == 1 && (env[parentIdx] == 1 && env[childIdx] == 2)
      # .mutate once and only once when front touches the hotspot
      mutID = 2
      willMutate = false
      mutationCount += 1
    end

    # .shift group affliation to match environment
    if mutID == 1
      groupID = env[childIdx] == 1 ? 1 : 2
    else
      groupID = env[childIdx] == 1 ? 3 : 4
    end

    # .fill new node, add to front if viable
    node = graph[childIdx]
    node.filled = true
    node.time = time
    node.ID_1 = groupID
    node.ID_2 = graph[parentIdx].ID_2
    node.ID_3 = mutID
    node.ancestor = parentIdx

    # .keep mutation number from parent
    if mutID == 2
      node.ID_4 = graph[parentIdx].ID_4 == 0 ? mutationCount : graph[parentIdx].ID_4
    end

    # .don't add site to front if it already has no empty nbors
    if node.nbors > 0
      add!(active[groupID], childIdx)
    end

    # .iterate through neighbors, updating nbor counts
    for nbor in view(cntns, :, childIdx)
      nbor == 0 && continue

      # .subtract one from all neighbors
      graph[nbor].nbors -= 1

      # .if this site, or neighbor, is surrounded, then remove it
      if graph[nbor].nbors == 0 && graph[nbor].filled
        remove!(active[graph[nbor].ID_1], nbor)
      end
    end

    # .stop recording data for later analysis
    if childIdx >= stopCondition && isempty(front)
      front = vcat(active[1].data, active[2].data, active[3].data, active[4].data)
      continueRecording = false
    end
  end
  return (front=front, time=time, extinction=tExtinction, width=width, height=height)
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

function setPath(cli)::String
  opts = [
    "env_type", "initial_type", "width", "height", "selection", "mutation", "intensity", "radius", "density", "rngSeed"
  ]
  opts = Symbol.(opts)

  parsedOpts = []
  for opt in opts
    haskey(cli, Symbol(opt)) || continue
    val = getfield(cli, Symbol(opt))
    if eltype(val) <: Float64
      val = round.(val, digits=3)
    end
    push!(parsedOpts, "$(opt)_$(val)")
  end

  # optional paramters
  cli.gap > 0 && push!(parsedOpts, "gap_$(cli.gap)")
  cli.env_type == "circle" && push!(parsedOpts, "sep_$(cli.separation)")
  cli.detailed_analytics && push!(parsedOpts, "da")
  cli.standing_variation && push!(parsedOpts, "sv")
  cli.animate && push!(parsedOpts, "animated")

  path = mkpath(cli.outputPath * "/" * join(parsedOpts, ","))

  if ispath(path) && ~isempty(readdir(path))
    if cli.rewrite
      foreach(rm, filter(endswith(".txt"), readdir(path; join=true)))
      foreach(rm, filter(endswith(".jld2"), readdir(path; join=true)))
      return path
    else
      @info " ! output path not empty"
      exit()
    end
  end
  return path
end

end
