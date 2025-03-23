using JLD2, FileIO, Arrow
import StatsBase.Random: seed!
include("../../src/base/dataModels.jl")
include("../../src/base/model.jl")
include("../../src/base/environment.jl")
include("../../src/base/trees.jl")
include("../../src/calculations/scaling.jl")

using CairoMakie, ColorSchemes
import ColorSchemes: tab10, ColorScheme
import Colors: RGB

using .Observables, .Model, .GenealogicalTree
using Parameters

module DataContainers
@kwdef mutable struct Results
  heatmap::Union{Missing,Vector{Float32}} = missing
  domainAreas::Union{Missing,Vector{Float32}} = missing
  htspts::Union{Missing,Vector{NTuple{2,Int64}}} = missing
  msd::Union{Missing,Vector{Float32}} = missing
  frontVariance::Union{Missing,Vector{Float32}} = missing
  sectorVariance::Union{Missing,Vector{Float32}} = missing
  mutationalFreq::Union{Missing,Vector{Float32}} = missing
  sectorSize::Union{Missing,Vector{Float32}} = missing
  sectorMV::Union{Missing,Matrix{Float32}} = missing
end
end

function images(graph, env, width, height, cli, bndry=missing)
  for label in [:ID_1, :ID_2, :ID_3, :ID_4]
    show_which == label || continue
    labelField = reshape(getfield.(graph, label), (width, height))

    fig = Figure(; size=(500, 0.87 * 500 * height / width))
    ax = Axis(fig[1, 1])
    hidedecorations!(ax)

    if (label == :ID_2 || label == :ID_4)
      colorrange = (1, width)
      cmap = tab10
    elseif label == :ID_1
      colorrange = (1, 6)
      cmap = ColorSchemes.gist_earth
    else
      colorrange = (1, 3)
      cmap = ColorScheme([RGB(1, 0, 0), RGB(0, 0, 0), RGB(1, 1, 0)])
    end

    #. draw proper ID's 
    heatmap!(ax, labelField; colormap=cmap, colorrange=colorrange, lowclip=:white, highclip=:yellow)

    #. draw hotspots if non zero intensity
    if cli.intensity > 0 && label != :ID_1
      heatmap!(ax, reshape(env, (width, height)); colorrange=(1.1, 1.2), lowclip=:transparent, highclip=(:blue, 0.4))
    end

    if label == :ID_3 && contains(cli.initial_type, "split") && cli.standing_variation
      scatter!(ax, bndry, collect(1:length(bndry)); color=:dodgerblue, markersize=4)
    end
    display(fig)
  end
end

function main!(graph, trial, dataModels, cli, env, results)
  # .reset models
  resetModels!(dataModels)

  # .reset graph properties
  resetGraph!(graph, cntns)

  snapshot = simulate!(graph, cntns, cli, dataModels, env)
  sources = Model.graphSources(cli)

  if cli.detailed_analytics
    # .determine the phylo of the periphery population 
    #phylo, branchPoints, branchTimes = GenealogicalTree.genealogy(graph, sources)
    #geneticLabels = getfield.(graph[sources], :ID_2)
    #pids = getfield.(graph[sources], :ID_2)
    #mutids = getfield.(graph[sources], :ID_3)

    #results.sectorSize = DataModels.lateralSectorSize(graph, opts.dims)

    #ng = JLD2.Group(file, "trial_$trial")
    #foreach(kv -> ng[String(first(kv))] = last(kv), DataModels.unpack(dataModels))
  end

  if contains(cli.initial_type, "split") && cli.standing_variation
    # .trace boundary and determine scaling
    bndry = Observables.traceBoundary(graph, cli.width, cli.height; window=20, delim=3)

    #. update the sector mean and variance
    if trial == 1
      results.sectorMV = zeros(Float32, length(bndry), 3)
      results.sectorMV[:, 1] .= 1
      results.sectorMV[:, 2] .= bndry 
    else
      Observables.onlineUpdates!(results.sectorMV, bndry)
    end

    #results.msd = Observables._lateralMSD(bndry)
    #results.sectorVariance = Observables._roughening(bndry)
    #results.mutationalFreq = Observables._mutantFrequency(graph, cli.width, cli.height)
  end

  # .save snapshop of simulation
  if trial == 1
    if @isdefined(bndry)
      images(graph, env, cli.width, cli.height, cli, bndry)
    else
      images(graph, env, cli.width, cli.height, cli, missing)
    end
  end

  # .domain size scaling
  ξ = (0.0, 0.0)
  if cli.mutation > 0
    ξ, domain_hist = Observables.domainSizes(graph, cli.dims; nbins=10000)
    if ismissing(results.domainAreas)
      results.domainAreas = domain_hist
    else
      length(domain_hist) == length(results.domainAreas) && (results.domainAreas .+= domain_hist)
    end
  end

  # .surviving numbers
  survFreq = Observables.survivalFreq(graph, sources, :ID_3; species=3)

  return round.((survFreq..., snapshot.time, snapshot.extinction, ξ...), digits=3)
end

# .command line options (input)
cli = (
  width=500,
  height=500,
  dims=(500, 1000),
  initial_type="split",
  standing_variation=true,
  mutation=0.01,
  detailed_analytics=true,
  animate=false,
  selection=0.1,
  compensation=0.1,
  rngSeed=4,
  landscape=nothing,
  env_type="uniform",
  intensity=0,
  density=0.0,
  radius=10,
  gap=0,
  numberTrials=100,
  numberSamples=10,
  singleMutant=false
)

# .set seed for reproducabilty
rngSeed = cli.rngSeed == 0 ? rand() : cli.rngSeed
seed!(rngSeed)

# .create containers for Observables and associated functions
dataModels = emptyModels(; animation=cli.animate)

# .create landscape of hotspots
env, htspts, cli = createEnvironment!(cli)

# .assert proper bounds of input parameters, c.g. Castillo 2020
if ~(0 <= cli.selection < 1 && 0 <= cli.compensation <= cli.selection)
  println("check input parameters\n")
  @warn "check input parameters\n"
end

# .send command line parameters to stdout and mkdirs
cli = merge(cli, (rngSeed=rngSeed,))

# .dimensions
lx = cli.width
ly = cli.height

# .construct graphs
periodic = contains(cli.initial_type, "split") ? false : true
graph, cntns = hexGraph((lx, ly), periodic)

results = DataContainers.Results()

show_which = :ID_3

for i in 1:cli.numberTrials
  data = main!(graph, i, dataModels, cli, env, results)
end

s_mean, _, s_var = Observables.retrieveOnlineUpdates(results.sectorMV) 

#. plot our results
fig, ax = scatter(collect(1:ly), s_mean, color=:blue, label="mean", markersize=2)
ylims!(ax, -10, 10)
display(fig)
fig, ax = scatter(collect(1:ly), sqrt.(s_var), color=:blue, label="roughness")