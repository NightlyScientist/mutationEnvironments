module Model
using StatsBase
include("../base/containers.jl")
import ArgParse: ArgParseSettings, parse_args, add_arg_table!, parse_item, add_arg_group!
import StaticArrays: SVector
import .HashMaps: HashVec, add!, remove!
export simulate!, resetGraph!, parseArgs, setPath, hexGraph

parse_item(::Type{NTuple{3,T}}, x::AbstractString) where {T} = Tuple(convert.(T, parse.(Float64, split(x, ','))))

function parseArgs(printinfo=true)
  addArgs!(sts, name; kwargs...) = add_arg_table!(sts, "--$(name)", Dict(kwargs))

  sts = ArgParseSettings()
  addArgs!(sts, "landscape"; required=false, arg_type=String)
  addArgs!(sts, "env_type"; arg_type=String, default="uniform")
  addArgs!(sts, "separation"; arg_type=Int)
  addArgs!(sts, "gap"; arg_type=Int, default=0)
  addArgs!(sts, "neighbors"; arg_type=Int, default=1)

  addArgs!(sts, "numberTrials"; arg_type=Int, default=1)
  addArgs!(sts, "numberSamples"; arg_type=Int, default=50)

  addArgs!(sts, "width"; required=true, arg_type=Int64)
  addArgs!(sts, "height"; required=true, arg_type=Int64)
  addArgs!(sts, "selection"; required=true, arg_type=Float64)
  addArgs!(sts, "compensation"; required=true, arg_type=Float64)
  addArgs!(sts, "mutation"; required=true, arg_type=Float64)
  addArgs!(sts, "intensity"; required=true, arg_type=Float64)
  addArgs!(sts, "density"; required=false, arg_type=Float64, default=0.09)
  addArgs!(sts, "radius"; required=false, arg_type=Int64, default=10)

  addArgs!(sts, "rngSeed"; arg_type=Int, default=1)
  addArgs!(sts, "outputPath"; required=true, arg_type=String)

  addArgs!(sts, "heatmap"; action=:store_true)
  addArgs!(sts, "animate"; action=:store_true)
  addArgs!(sts, "rewrite"; action=:store_true)
  addArgs!(sts, "printInfo"; action=:store_true)
  addArgs!(sts, "detailed_analytics"; action=:store_true)
  addArgs!(sts, "initial_type"; arg_type=String, default="alt")

  # exclusive option groups
  add_arg_group!(sts, "exlusive"; exclusive=true)
  addArgs!(sts, "standing_variation"; action=:store_true)
  addArgs!(sts, "singleMutant"; action=:store_true)
  parsedArgs = namedtuple(parse_args(sts))

  if printinfo
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

mutable struct Status
  mutationCount::UInt64
  extinction_time::Float64
  has_mutated::Bool
  only_one_mutation::Bool
  hard_obstacles::Bool
  terminal_index::Int64
end

graphSources(cli) = collect(1:(cli.width)) .+ cli.width * (cli.height - 1)

function populate!(graph, dims, env, opts; row=1, num=3, standingVar=false)
  if standingVar
    if contains(opts.initial_type, "split")
      # right half is population 3
      strainID = ones(Int, dims[1])
      strainID[cld(dims[1], 2) + 1:end] .= 3
    elseif contains(opts.initial_type, "alt")
      # alt wild-type and population 3
      types = [1, 3]
      strainID = [types[mod(x, 1:2)] for x in 1:dims[1]]
    else
      # random initial population
      strainID = rand([1, 3], dims[1])
    end
  else
    # all wild-type, to be used with mutation rates
    strainID = ones(Int, dims[1])
  end

  for col in 1:dims[1]
    #. set group affliliation with strain type (f, s, b), matching the environment
    nodeID = strainID[col]

    #. shift group affliation to match environment 
    env[col] == 2 && (nodeID += 3)

    nodeIndx = dims[1] * (row - 1) + col

    graph[nodeIndx].filled = true
    graph[nodeIndx].ID_1 = nodeID
    graph[nodeIndx].ID_2 = col
    graph[nodeIndx].ID_3 = strainID[col]
  end
  #@info "population initialized"
  #@info getfield.(graph[:, 1], :ID_1)
end

function hexGraph(dims)
  lx, ly = dims
  graph = Vector{Node}(undef, lx * ly)
  foreach(i -> graph[i] = Node(6), eachindex(graph))
  return graph
end

function ssa(rates)::Tuple{Int,Float64}
  acculm = accumulate(+, rates)
  η = rand() * last(acculm)
  for i in eachindex(acculm)
    η <= acculm[i] && return (i, rates[i])
  end
end

@inline function pbc(lx::Int64, x::Int64)
  if x < 1
    x = lx + x
  elseif x > lx
    x = x - lx
  end
  return x
end

function hardWall(lx::Int64, x::Int64)
  if x < 1
    return 1
  elseif x > lx
    return lx
  end
  return x
end

function newSite(rates, layer_previous, neighbors::Int64, index::Int64, periodic::Bool)
  lx = size(layer_previous, 1)
  if periodic
    n_neighbors = 2 * neighbors + 1
    _range = (-neighbors):neighbors
  else
    xlow = hardWall(lx, index - neighbors)
    xhigh = hardWall(lx, index + neighbors)
    _range = (xlow - index):(xhigh-index)
    #println("xlow: $xlow, xhigh: $xhigh")
    n_neighbors = xhigh - xlow + 1
  end
  subset_g = zeros(Float64, n_neighbors)
  subset = zeros(Int64, n_neighbors)

  counter = 0
  for dx in _range
    counter += 1
    x = pbc(lx, index + dx)
    ID_1 = layer_previous[x].ID_1
    subset_g[counter] = ID_1 > 0 ? rates[ID_1] : 0
    subset[counter] = x
    #print("$ID_1,")
  end
  #println("\n")
  # check that not all cells have are in a hard obstacle
  if all(subset_g .== 0)
    return -1
  end

  index, R = ssa(subset_g)
  return subset[index]
end

function createGeneration!(graph, rates, mut, env, t::Int64, neighbors::Int64, status::Status, periodic::Bool=false)
  layer_current = view(graph, :, t)
  layer_previous = view(graph, :, t - 1)

  #println("previous layer: $(t-1)")
  #println(getfield.(layer_previous, :ID_1))
  #@info getfield.(layer_previous, :ID_1)

  for index in eachindex(layer_current)

    # fill new site at index by sampling previous generation
    #@info "replicating now into $index"
    replicating_index = newSite(rates, layer_previous, neighbors, index, periodic)

    # don't fill if hard obstacle
    if replicating_index < 0
      #@info "replicating index: $replicating_index into $index)"
      layer_current[index].filled = false
      layer_current[index].ID_1 = 0
      continue
    end

    replicating_node = layer_previous[replicating_index]
    #@info "replicating index: $replicating_index into $index using $(layer_previous[replicating_index].ID_1)"

    # mutate childID based on mutation rate
    strainID = layer_previous[replicating_index].ID_3
    if !status.only_one_mutation && strainID == 1 && rand() < mut
      strainID = 2 # mutate in to mutant 1 -> 2
      status.mutationCount += 1
    end

    parent_env_type = replicating_node.ID_1
    child_env_type = env[index, t]
    if status.only_one_mutation && !status.has_mutated && strainID == 1 && (parent_env_type == 1 && child_env_type == 2)
      # .mutate once and only once when front touches the hotspot
      strainID = 2
      status.has_mutated = true
      status.mutationCount += 1
    end

    # inhereit group affliliation from ancestor
    # groupID = graph[parentIdx].ID_1
    #groupID = strainID

    # .shift group affliation to match environment
    groupID = strainID
    child_env_type == 2 && (groupID = strainID + 3)

    # fill new node
    layer_current[index].filled = true
    layer_current[index].ID_1 = groupID
    layer_current[index].ID_2 = replicating_node.ID_2
    layer_current[index].ID_3 = strainID
    layer_current[index].ancestor = replicating_index
    layer_current[index].time = t

    # keep mutation number from parent
    if strainID == 2
      layer_current[index].ID_4 = replicating_node.ID_4 == 0 ? status.mutationCount : replicating_node.ID_4
      #node.ID_4 = graph[parentIdx].ID_4 == 0 ? mutationCount : graph[parentIdx].ID_4
    end
  end
end

function populationCounts(graph, t)
  layer = view(graph, :, t)
  counts = zeros(Int, 6)

  for index in eachindex(layer)
    ID = layer[index].ID_1
    #print("$ID,")
    ID == 0 && continue
    counts[ID] += 1
  end
  #println(" ")
  return counts
end

function simulate!(graph, cli, dataModels, env; periodic=true)
  sel, ν, comp, mut = cli.selection, cli.intensity + 1, cli.compensation, cli.mutation
  width, height = cli.width, cli.height
  nbors = cli.neighbors

  #. wild type | mutant | bystander | wild type (env) | mutant (env) | bystander (env)
  rates = SVector{6,Float64}([1, 1 - sel, 1 - sel + comp, ν, ν * (1 - sel), ν * (1 - sel + comp)])

  #. check if env features are hard_obstacles (no growth)
  hard_obstacles = ν <= 0

  # .track number of mutations that occur
  singleMutant = cli.singleMutant || cli.env_type == "circle"
  status = Status(0, -1.0, false, singleMutant, hard_obstacles, cli.height)

  populate!(graph, (width, height), env, cli; num=6, standingVar=cli.standing_variation)

  # reshape the graph object to a 2D array
  graph = reshape(graph, (width, height))
  env = reshape(env, (width, height))

  #for i in 1:width
  #  print(graph[i, 1].ID_1)
  #  print(", ")
  #end

  #println("starting the main thing")
  populationCounts(graph, 1)

  @inbounds for t in 2:(cli.height)
    # update the graph with the next generation
    createGeneration!(graph, rates, mut, env, t, nbors, status, periodic)

    # fetch the population counts
    population_counts = populationCounts(graph, t)

    #. terminate simulation when bystander population is extinct
    if population_counts[3] + population_counts[6] == 0
      status.extinction_time < 0 && (status.extinction_time = t)
      status.terminal_index = t
      cli.heatmap || break
    end
  end

  front = graph[:, status.terminal_index]

  graph = reshape(graph, width * height)
  env = reshape(env, width * height)
  return (
    front=front,
    time=cli.height,
    extinction=status.extinction_time,
    width=width,
    height=height,
    stop_index=status.terminal_index
  )
end

resetGraph!(graph) = @simd for i in eachindex(graph)
  #resetGraph!(graph, cntns) = @inbounds @simd for i in eachindex(graph)
  graph[i].filled = false
  graph[i].ancestor = 0
  graph[i].ID_1 = 0
  graph[i].ID_2 = 0
  graph[i].ID_3 = 0
  graph[i].ID_4 = 0
  graph[i].time = 0.0
  graph[i].nbors = 0
end

function setPath(cli)::String
  opts = [
    "env_type",
    "initial_type",
    "width",
    "height",
    "neighbors",
    "selection",
    "compensation",
    "mutation",
    "intensity",
    "radius",
    "density",
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

  path = mkpath(cli.outputPath * "/1D/" * join(parsedOpts, ","))

  #. write ARGS to log file
  log_dir = joinpath(path, "logs")
  isdir(log_dir) || mkdir(log_dir)

  log_file = joinpath(log_dir, "log.md")
  open(log_file, "w") do file
    write(file, "command line: $(join(ARGS, ' '))\n\n")
    write(file, "Simulation Info:\n")
    foreach(k -> write(file, "  > $k  =>  $(cli[k])\n"), keys(cli))
  end

  if ispath(path) && ~isempty(readdir(path))
    if cli.rewrite
      foreach(rm, filter(endswith(".png"), readdir(path; join=true)))
      foreach(rm, filter(endswith(".jld2"), readdir(path; join=true)))
      foreach(rm, filter(endswith(".arrow"), readdir(path; join=true)))
      foreach(rm, filter(endswith(".csv"), readdir(path; join=true)))
      return path
    else
      @info " ! output path not empty"
      exit()
    end
  end
  return path
end

end
