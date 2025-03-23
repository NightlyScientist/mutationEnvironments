module Model
include("../../base/containers.jl")
include("../../common/settings.jl")
include("../../common/indexTools.jl")
include("../../base/grids.jl")

using Statistics
import StaticArrays: SVector
import .HashMaps: HashVec, add!, remove!
import .SimulationSettings: Setting, Settings, parse_settings, get_abbrs, print_settings, ensurePath
export simulate!, resetGraph!, hexGraph

function settings_check(parsedArgs::NamedTuple)
  parsedArgs.n_species <= 3 || error("n_species must be less than or equal to 3")
  parsedArgs.printInfo && print_settings(parsedArgs)
end

function simulation_settings()
  sts = Settings(; program_name="kpz_with_renewal")
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
    Setting(; name="initial_type", abbr="IT", arg_type=String, default="alt"),
    Setting(; name="fraction_remove", abbr="FR", arg_type=Float64, default=0.1),
    Setting(; name="frames", abbr="F", arg_type=Int64, default=10)
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
    # remove n % of the sites and add them to the active list
    rand() <= opts.fraction_remove && continue

    x, y = reconstructCoordinates(i, dims...)
    band_width = 10
    nodeID = floor(x / band_width) % 2 == 0 ? 1 : 2
    _nodeID = floor(x / band_width) % 2 == 0 ? 1 : 2

    #if floor(y / band_width) % 2 == 0
    #  nodeID = floor(x / band_width) % 2 == 0 ? 2 : 1
    #  _nodeID = floor(x / band_width) % 2 == 0 ? 2 : 1
    #end

    if opts.mutation > 0.0 || ~opts.standing_variation
      nodeID = 1
      _nodeID = 1
    end

    # shift group affliation to match environment
    env[i] == 2 && (nodeID += n_species)

    nodeIndx = i
    graph[nodeIndx].filled = true
    graph[nodeIndx].ID_1 = nodeID
    graph[nodeIndx].ID_2 = i % dims[1]
    #graph[nodeIndx].ID_3 = floor(x / band_width) % 2 == 0 ? 1 : 2
    graph[nodeIndx].ID_3 = _nodeID
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
  r = zeros(Float64, length(active))
  for i in eachindex(active)
    r[i] = length(active[i].data) * rates[i]
  end

  acculm = accumulate(+, r)
  η = rand() * last(acculm)
  for i in eachindex(acculm)
    η <= acculm[i] && return (i, rates[i])
  end
end

"""remove x-fraction of the population"""
function remodel_landscape!(active, graph, cntns, population_counts, cli, time_scale)
  for i in eachindex(graph)
    graph[i].filled || continue
    graph[i].nbors == 0 || continue

    # rand() <= cli.fraction_remove || continue
    rand() <= (1 - exp(-(time - graph[i].time) / time_scale)) || continue

    # remove node from active list
    removed = remove!(active[graph[i].ID_1], i)
    removed && (population_counts[graph[i].ID_1] -= 1)

    graph[i].filled = false
    graph[i].ancestor = 0
    graph[i].ID_1 = 0
    graph[i].ID_2 = 0
    graph[i].ID_3 = 0
    graph[i].ID_4 = 0
    graph[i].time = 0.0
    graph[i].nbors = 0

    # insert neighbors around the site into the active list
    for nbor in view(cntns, :, i)
      nbor == 0 && continue
      graph[nbor].nbors += 1

      if graph[nbor].filled
        added = add!(active[graph[nbor].ID_1], nbor)
        added && (population_counts[graph[nbor].ID_1] += 1)
      else
        graph[i].nbors += 1
      end
    end
  end
end

mutable struct State
  mutationCount::Int
  tExtinction::Float64
  willMutate::Bool
  singleMutant::Bool
  hard_obstacles::Bool
  mut::Float64
  time::Float64
  stop_index::Int64
  num_pop::Int64
  rates::SVector
end

function setup_system_state(cli)::State
  sel, intensity, mut, comp = cli.selection, cli.intensity, cli.mutation, cli.compensation
  ν = cli.intensity + 1

  # check if env features are hard_obstacles (no growth)
  hard_obstacles = ν <= 0

  # wild type | mutant | bystander | wild type (env) | mutant (env) | bystander (env)
  rates = SVector{6,Float64}([1, 1 - sel, 1 - sel + comp, ν, ν * (1 - sel), ν * (1 - sel + comp)])

  if cli.n_species == 2
    # wild type | mutant | wild type (env) | mutant (env)
    rates = SVector{4,Float64}([1, 1 + intensity, 1 - sel, (1 - sel) * (1 + intensity)])
  end

  singleMutant = cli.singleMutant || cli.env_type == "circle"
  time = 0.0
  num = cli.n_species * 2

  # .save snapshot of the front as it touches the top
  width, height = cli.width, cli.height
  stop_index = width * height
  return State(0, -1, true, singleMutant, hard_obstacles, mut, 0.0, stop_index, num, rates)
end

function replicate_active!(active, rates, env, graph, cntns, population_counts, state, cli)
  number_active = sum(population_counts)

  for i in 1:number_active
    # all sites should be filled now
    sum(population_counts) == 0 && break

    # find which cell will grow next
    groupID, R = ssa(active, rates)

    # parent and child index
    isempty(active[groupID].data) && (@info "empty active list"; break)
    parentIdx = rand(active[groupID].data)
    childIdx = nextIndex(graph, cntns, parentIdx)

    # update time
    state.time += -log(rand()) / R

    # .mutate childID based on mutation rate, and increment mutations counts
    strainID = graph[parentIdx].ID_3
    !state.singleMutant && strainID == 1 && rand() < state.mut && (strainID = 2; state.mutationCount += 1)

    if state.singleMutant && state.willMutate && strainID == 1 && (env[parentIdx] == 1 && env[childIdx] == 2)
      # .mutate once and only once when front touches the hotspot
      strainID = 2
      state.willMutate = false
      state.mutationCount += 1
    end

    # .shift group affliation to match environment
    env[childIdx] == 2 ? (groupID = strainID + cli.n_species) : (groupID = strainID)

    # .fill new node, add to front if viable
    node = graph[childIdx]
    node.filled = true
    node.time = state.time
    node.ID_1 = groupID
    node.ID_2 = graph[parentIdx].ID_2
    node.ID_3 = strainID
    node.ancestor = parentIdx

    # .keep mutation number from parent
    if strainID == 2
      node.ID_4 = graph[parentIdx].ID_4 == 0 ? state.mutationCount : graph[parentIdx].ID_4
    end

    #. check if child is in a hard obstacle
    good_site = (env[childIdx] == 2 && state.hard_obstacles) ? false : true

    # don't add site to front if it already has no empty nbors if node.nbors > 0
    if node.nbors > 0 && good_site
      add!(active[groupID], childIdx)
      population_counts[groupID] += 1
    end

    # .iterate through neighbors, updating nbor counts
    for nbor in view(cntns, :, childIdx)
      nbor == 0 && continue

      # .subtract one from all neighbors
      graph[nbor].nbors -= 1

      # .if this site, or neighbor, is surrounded, then remove it
      if graph[nbor].nbors == 0 && graph[nbor].filled
        removed = remove!(active[graph[nbor].ID_1], nbor)
        removed && (population_counts[graph[nbor].ID_1] -= 1)
      end
    end
  end
end

function simulate!(graph, cntns, cli, dataModels, env)
  state::State = setup_system_state(cli)

  population_counts = zeros(Int, state.num_pop)
  #active = populate!(
  #  graph, cntns, (width, height), env, cli, population_counts; num=num, standingVar=cli.standing_variation
  #)

  dims = (state.width, state.height)
  active = populate_fill!(
    graph, cntns, dims, env, cli, population_counts; num=state.num_pop, standingVar=cli.standing_variation
  )

  time_scale = 1 / sum(state.rates) * 100

  @inbounds for _ in 1:(cli.frames)
    #sum(population_counts) == 0 && (@info "population counts is now zero"; break)
    #sum(population_counts[1:(cli.n_species)]) == 0 && tExtinction < 0 && (tExtinction = time)

    # replicate all active sites
    replicate_active!(active, state.rates, env, graph, cntns, population_counts, state, cli)

    # record data for each frame
    for (_, measure) in dataModels
      push!(measure.data, measure.func(active, graph, time))
    end

    # remodel the landscape by removing half of the filled sites
    remodel_landscape!(active, graph, cntns, population_counts, cli, time_scale)
  end

  front = reduce(vcat, getfield.(active, :data))
  return (front=front, time=time, extinction=state.tExtinction, stop_index=state.stop_index)
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
end
