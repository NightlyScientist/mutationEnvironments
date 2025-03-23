module Model
using Statistics
include("../base/containers.jl")
include("../common/settings.jl")
include("../base/grids.jl")
include("initializations.jl")
import StaticArrays: SVector
import .HashMaps: HashVec, add!, remove!
import .SimulationSettings: Setting, Settings, parse_settings, get_abbrs, print_settings, ensurePath

function settings_check(parsedArgs::NamedTuple)
  parsedArgs.n_species <= 3 || error("n_species must be less than or equal to 3")
  parsedArgs.printInfo && print_settings(parsedArgs)
end

function simulation_settings()
  sts = Settings()
  sts.exclusive = [("standing_variation", "sv"), ("singleMutant", "sm")]
  sts.flags = [("heatmap", "hm"), ("animate", "ani"), ("rewrite", ""), ("printInfo", ""), ("detailed_analytics", "da"), ("fixed_initializations", "fxinit")]
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
  @inbounds for i in eachindex(active)
    r[i] = length(active[i].data) * rates[i]
  end

  acculm = accumulate(+, r)
  η = rand() * last(acculm)
  for i in eachindex(acculm)
    η <= acculm[i] && return (i, rates[i])
  end
end

# mutable struct SimulationState
# end

function simulate!(graph, cntns, cli, dataModels, env, save_state)
  (; width, height) = cli
  (sel, intensity, mut, comp) = cli.selection, cli.intensity, cli.mutation, cli.compensation

  #. check if env features are hard_obstacles (no growth)
  ν = cli.intensity + 1
  hard_obstacles = ν <= 0

  if cli.n_species == 2
    #. wild type | mutant | wild type (env) | mutant (env)
    rates = SVector{4,Float64}([1, 1 - sel, ν, (1 - sel) * ν])
  else
    #. wild type | mutant | bystander | wild type (env) | mutant (env) | bystander (env)
    rates = SVector{6,Float64}([1, 1 - sel, 1 - sel + comp, ν, ν * (1 - sel), ν * (1 - sel + comp)])
  end

  time = 0.0

  population_counts = zeros(Int, cli.n_species * 2)
  active = initializePopulation!(graph, cntns, env, cli, population_counts, save_state)

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

  @inbounds while true
    sum(population_counts) == 0 && break
    sum(population_counts[1:(cli.n_species)]) == 0 && tExtinction < 0 && (tExtinction = time)

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
    strainID = graph[parentIdx].ID_3
    if !singleMutant && strainID == 1 && rand() < mut
      strainID = 2
      mutationCount += 1
    end

    if singleMutant && willMutate && strainID == 1 && (env[parentIdx] == 1 && env[childIdx] == 2)
      # .mutate once and only once when front touches the hotspot
      strainID = 2
      willMutate = false
      mutationCount += 1
    end

    # .shift group affliation to match environment
    groupID = env[childIdx] == 2 ? strainID + cli.n_species : strainID

    # .fill new node, add to front if viable
    node = graph[childIdx]
    node.filled = true
    node.time = time
    node.ID_1 = groupID
    node.ID_2 = graph[parentIdx].ID_2
    node.ID_3 = strainID
    node.ancestor = parentIdx

    # .keep mutation number from parent
    if strainID == 2
      node.ID_4 = graph[parentIdx].ID_4 == 0 ? mutationCount : graph[parentIdx].ID_4
    end

    #. check if child is in a hard obstacle
    good_site = (env[childIdx] == 2 && hard_obstacles) ? false : true

    # .don't add site to front if it already has no empty nbors if node.nbors > 0
    if node.nbors > 0 && good_site
      add!(active[groupID], childIdx)
      population_counts[groupID] += 1
    end

    # .iterate through neighbors, updating nbor counts
    updateNeighbors!(graph, cntns, childIdx, active, population_counts)

    # .stop recording data for later analysis
    if childIdx >= stopCondition && isempty(front)
      front = reduce(vcat, getfield.(active, :data))
      continueRecording = false
      #cli.heatmap || break
    end

    #.terminate simulation when bystander population is extinct
    if population_counts[cli.n_species] + population_counts[cli.n_species * 2] == 0
      #if (length(active[3].data) + length(active[6].data)) == 0
      if tExtinction < 0
        tExtinction = time
        stop_index = childIdx
        front = reduce(vcat, getfield.(active, :data))
      end
      cli.heatmap || break
    end
  end
  return (front=front, time=time, extinction=tExtinction, stop_index=stop_index)
end
end
