module Model
using Statistics
include("../../base/containers.jl")
include("../../common/settings.jl")
include("../../common/indexTools.jl")
include("../../base/grids.jl")

import StaticArrays: SVector
import .HashMaps: HashVec, add!, remove!
import .SimulationSettings: Setting, Settings, parse_settings, get_abbrs, print_settings, ensurePath
export simulate!, resetGraph!, hexGraph

function settings_check(parsedArgs::NamedTuple)
  parsedArgs.n_species <= 3 || error("n_species must be less than or equal to 3")
  if parsedArgs.printInfo
    println("Simulation Info:")
    foreach(k -> println("  > $k  =>  $(parsedArgs[k])"), keys(parsedArgs))
  end
end

function simulation_settings()
  sts = Settings(; program_name = "with_replacement")
  sts.exclusive = [("standing_variation", "sv"), ("singleMutant", "sm")]
  sts.flags = [("heatmap", "hm"), ("animate", "ani"), ("rewrite", ""), ("printInfo", ""), ("detailed_analytics", "da")]
  sts.options = [
    Setting(; name="landscape", arg_type=String),
    Setting(; name="env_type", abbr="ENV", arg_type=String, default="uniform"),
    Setting(; name="separation", abbr="Sep", arg_type=Int),
    Setting(; name="gap", abbr="gap", arg_type=Int, default=0),
    Setting(; name="numberTrials", arg_type=Int, default=1),
    Setting(; name="numberSamples", arg_type=Int, default=50),
    Setting(; name="width", abbr="LX", required=true, arg_type=Int64),
    Setting(; name="height", abbr="LY", required=true, arg_type=Int64),
    Setting(; name="compensation", abbr="C", required=false, arg_type=Float64, default=0.0),
    Setting(; name="selection", abbr="S", required=false, arg_type=Float64, default=0.0),
    Setting(; name="mutation", abbr="M", required=false, arg_type=Float64, default=0.0),
    Setting(; name="intensity", abbr="I", required=false, arg_type=Float64, default=0.0),
    Setting(; name="density", abbr="D", required=false, arg_type=Float64, default=0.00),
    Setting(; name="radius", abbr="R", required=false, arg_type=Int64, default=10),
    Setting(; name="n_species", abbr="N", arg_type=Int64, default=2),
    Setting(; name="rngSeed", abbr="RS", arg_type=Int, default=1),
    Setting(; name="outputPath", required=true, arg_type=String),
    Setting(; name="initial_type", abbr="IT", arg_type=String, default="alt")
  ]
  parsed = parse_settings(sts)
  parsed = ensurePath(sts, parsed)
  settings_check(parsed)
  return parsed
end

graphFront(cli) = collect(1:(cli.width)) .+ cli.width * (cli.height - 1)

"""fill the entire landscape with a wild type and mutant"""
function populate_fill!(graph, cntns, dims, env, opts, container; row=1, num=3, standingVar=false)
  active = Vector{HashVec{UInt32}}(undef, num)
  foreach(i -> active[i] = HashVec{UInt32}(), collect(1:1:num))

  n_species = (num == 6) ? 3 : 2

  for i in eachindex(graph)
    # remove 20 % of the sites and add them to the active list
    #rand() <= 0.2 && continue

    x, _ = reconstructCoordinates(i, dims...)
    nodeID = floor(x / 10) % 2 == 0 ? 1 : 2

    # shift group affliation to match environment
    env[i] == 2 && (nodeID += n_species)

    nodeIndx = i
    graph[nodeIndx].filled = true
    graph[nodeIndx].ID_1 = nodeID
    graph[nodeIndx].ID_2 = i % dims[1]
    graph[nodeIndx].ID_3 = floor(x / 10) % 2 == 0 ? 1 : 2
    graph[nodeIndx].time = rand() * 100

    # don't add site to front if it already has no empty nbors
    if graph[nodeIndx].nbors > 0
      add!(active[nodeID], nodeIndx)
      container[nodeID] += 1
    end

    for nbor in view(cntns, :, nodeIndx)
      nbor == 0 && continue

      # subtract one from all neighbors
      graph[nbor].nbors -= 1

      # if this site, or neighbor, is surrounded, then remove it
      if graph[nbor].nbors == 0 && graph[nbor].filled
        remove!(active[graph[nbor].ID_1], nbor)
        container[graph[nbor].ID_1] -= 1
      end
    end
  end
  return active
end

function populate!(graph, cntns, dims, env, opts, container; row=1, num=3, standingVar=false)
  active = Vector{HashVec{UInt32}}(undef, num)
  foreach(i -> active[i] = HashVec{UInt32}(), collect(1:1:num))

  n_species = (num == 6) ? 3 : 2

  if standingVar
    if contains(opts.initial_type, "split")
      strainID = ones(Int, dims[1])
      strainID[(cld(dims[1], 2) + 1):end] .= n_species
    elseif contains(opts.initial_type, "alt")
      types = [1, n_species]
      strainID = [types[mod(x, 1:2)] for x in 1:dims[1]]
      #strainID = mod.(collect(1:dims[1]), 2) .+ 1
    else
      strainID = rand([1, n_species], dims[1])
    end
  else
    strainID = ones(Int, dims[1])
  end

  for col in 1:dims[1]
    #. set group affliliation with strain type (f, s, b), matching the environment
    nodeID = strainID[col]

    # shift group affliation to match environment
    env[col] == 2 && (nodeID += n_species)

    nodeIndx = dims[1] * (row - 1) + col
    graph[nodeIndx].filled = true
    graph[nodeIndx].ID_1 = nodeID
    graph[nodeIndx].ID_2 = col
    graph[nodeIndx].ID_3 = strainID[col]

    # don't add site to front if it already has no empty nbors
    if graph[nodeIndx].nbors > 0
      add!(active[nodeID], nodeIndx)
      container[nodeID] += 1
    end

    for nbor in view(cntns, :, nodeIndx)
      nbor == 0 && continue

      # subtract one from all neighbors
      graph[nbor].nbors -= 1

      # if this site, or neighbor, is surrounded, then remove it
      if graph[nbor].nbors == 0 && graph[nbor].filled
        remove!(active[graph[nbor].ID_1], nbor)
        container[graph[nbor].ID_1] -= 1
      end
    end
  end
  return active
end

function ssa(rates)::Tuple{Int,Float64}
  acculm = accumulate(+, rates)
  η = rand() * last(acculm)
  for i in eachindex(acculm)
    η <= acculm[i] && return (i, rates[i])
  end
end

"""iterate through the previouse time step's site's neighbors; return rates"""
function previousStepNeighbors(graph, cntns, ref_index, rates)
  filled_nbors = 0
  counter = 1
  nbors = zeros(Int, 7)
  _rates = zeros(Float64, 7)

  # neighbors
  for nbor in view(cntns, :, ref_index)
    nbor == 0 && continue
    graph[nbor].filled || continue

    filled_nbors += 1
    _rates[counter] = rates[graph[nbor].ID_1]
    nbors[counter] = nbor
    counter += 1
  end

  # self
  if graph[ref_index].filled
    _rates[7] = rates[graph[ref_index].ID_1]
    nbors[7] = ref_index
    filled_nbors += 1
  end

  filled_nbors == 0 && return -1
  return nbors[ssa(_rates)[1]]
end

function timestep!(
  max_index, current, previous, cntns, cli, env, mutationCount, willMutate, singleMutant, rates, mut, hard_obstacles
)
  #for site_index in 1:min(max_index, length(current))
  for site_index in eachindex(current)
    # get the site information from the previous time
    selected_index = previousStepNeighbors(previous, cntns, site_index, rates)

    selected_index == -1 && continue

    # .mutate childID based on mutation rate, and increment mutations counts
    strainID = previous[selected_index].ID_3
    if !singleMutant && strainID == 1 && rand() < mut
      strainID = 2
      mutationCount += 1
    end

    if singleMutant && willMutate && strainID == 1 && (env[selected_index] == 1 && env[site_index] == 2)
      # .mutate once and only once when front touches the hotspot
      strainID = 2
      willMutate = false
      mutationCount += 1
    end

    # .shift group affliation to match environment
    if env[site_index] == 2
      groupID = strainID + cli.n_species
    else
      groupID = strainID
    end

    # .fill new node, add to front if viable
    node = current[site_index]
    node.filled = true
    node.time = max_index
    node.ID_1 = groupID
    #node.ID_2 = rand(1:1000)
    node.ID_2 = previous[selected_index].ID_2
    node.ID_3 = strainID
    node.ancestor = selected_index

    # .keep mutation number from parent
    if strainID == 2
      node.ID_4 = previous[selected_index].ID_4 == 0 ? mutationCount : previous[selected_index].ID_4
    end

    #. check if child is in a hard obstacle, change site to unfilled
    good_site = (env[site_index] == 2 && hard_obstacles) ? false : true
    if !good_site
      node.filled = false
    end
  end
end

function simulate!(graph, cntns, cli, dataModels, env)
  sel = cli.selection
  intensity = cli.intensity
  mut = cli.mutation

  ν = cli.intensity + 1
  comp = cli.compensation

  #. check if env features are hard_obstacles (no growth)
  hard_obstacles = ν <= 0

  if cli.n_species == 2
    #. wild type | mutant | wild type (env) | mutant (env)
    rates = SVector{4,Float64}([1, 1 + intensity, 1 - sel, (1 - sel) * (1 + intensity)])
  else
    #. wild type | mutant | bystander | wild type (env) | mutant (env) | bystander (env)
    rates = SVector{6,Float64}([1, 1 - sel, 1 - sel + comp, ν, ν * (1 - sel), ν * (1 - sel + comp)])
  end

  width = cli.width
  height = cli.height

  time = 0.0

  num = cli.n_species * 2

  population_counts = zeros(Int, num)
  active = populate!(
    graph, cntns, (width, height), env, cli, population_counts; num=num, standingVar=cli.standing_variation
  )

  nrecord::Int64 = cld(width * height, cli.numberSamples)
  itrCntr::Int64 = nrecord
  continueRecording::Bool = true

  # .track number of mutations that occur
  mutationCount = 0

  # .save snapshot of the front as it touches the top
  front = Int32[]
  stopCondition = width * (height - 1) + 1
  stop_index = width * height

  tExtinction = -1

  willMutate = true
  singleMutant = cli.singleMutant || cli.env_type == "circle"

  # this is memory expensive, but should save time from allocations later
  alt_graph = deepcopy(graph)

  @inbounds for t in 1:(2 * height)
    if sum(population_counts[1:(cli.n_species)]) == 0 && tExtinction < 0
      tExtinction = time
    end

    # convert time into a y-index for searching
    next_row = width * (t + 1)

    # decide which container to use in updating, and what is considered 'previous'
    if t % 2 == 0 # use graph
      timestep!(
        next_row, graph, alt_graph, cntns, cli, env, mutationCount, willMutate, singleMutant, rates, mut, hard_obstacles
      )
    else #use alt_graph
      timestep!(
        next_row, alt_graph, graph, cntns, cli, env, mutationCount, willMutate, singleMutant, rates, mut, hard_obstacles
      )
    end
  end
  return (front=front, time=time, extinction=tExtinction, width=width, height=height, stop_index=stop_index)
end
end
