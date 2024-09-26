using JLD2, FileIO, Arrow
import StatsBase.Random: seed!
include("../calculations/scaling.jl")
include("../calculations/hotspotInfluence.jl")
include("dataModels.jl")
include("model.jl")
include("environment.jl")
include("trees.jl")

import CairoMakie: Figure, Axis, save, heatmap, heatmap!, hlines!, scatter!, hidedecorations!
import ColorSchemes: tab10, ColorScheme
import Colors: RGB

using .Observables, .Model, .GenealogicalTree

module DataContainers
mutable struct Results
  heatmap::Vector{Float32}
  domainAreas::Union{Missing,Vector{Float32}}
  influentialHotspots::Vector{UInt32}
  htspts::Vector{NTuple{2,Int64}}
end
end

function images(graph, width, height, cli, savePath, bndry=missing)
  for label in [:ID_1, :ID_2, :ID_3, :ID_4]
    colorrange = (label == :ID_2 || label == :ID_4) ? (1, width) : (1, 6)
    if (label == :ID_2 || label == :ID_4)
      colorrange = (1, width)
      cmap = tab10
    else
      colorrange = (1, 2)
      # ColorScheme([RGB(1, 0, 0), RGB(0, 0.5, 0.5), RGB(0, 1, 0), RGB(1, 1, 0), RGB(0, 1, 1), RGB(0, 0, 0)])
      cmap = ColorScheme([RGB(1, 0, 0), RGB(1, 1, 0)])
    end

    labelField = reshape(getfield.(graph, label), (width, height))
    if label == :ID_1
      # .don't draw hotspots if intensity is zero
      remapping = cli.intensity == 0.0 ? Int64[1, 1, 2, 2] : Int64[1, 3, 2, 3]
      labelField = map(t -> remapping[t], labelField)
    end
    fig = Figure(; size=(500, 0.87 * 500 * height / width))
    ax = Axis(fig[1, 1])
    heatmap!(ax, labelField; colormap=cmap, colorrange=colorrange, lowclip=:white, highclip=:black)

    if label == :ID_3 && contains(cli.initial_type, "split") && cli.standing_variation && cli.mutation == 0
      scatter!(ax, bndry, collect(1:length(bndry)); color=:dodgerblue, markersize=4)
    end
    hidedecorations!(ax)
    save(joinpath(savePath, "snapshots_$(label).png"), fig)
  end
end

function main!(graph, trial, dataModels, cli, env, _resutls, save_path)
  # .reset models
  resetModels!(dataModels)

  # .reset graph properties
  resetGraph!(graph, cntns)

  snapshot = simulate!(graph, cntns, cli, dataModels, env)
  sources = Model.graphSources(cli)

  if cli.detailed_analytics
    # .determine the phylo of the periphery population 
    phylo, branchPoints, branchTimes = GenealogicalTree.genealogy(graph, sources)
    geneticLabels = getfield.(graph[sources], :ID_2)

    jldopen(save_path * "/data_phylo.jld2", "a+") do file
      ng = JLD2.Group(file, "trial_$trial")
      ng["end_time"] = snapshot.time
      ng["source_labels"] = geneticLabels
      ng["phylogeny"] = phylo
      ng["branch_points"] = branchPoints
      ng["branch_times"] = branchTimes
    end

    # save dataModels
    jldopen(savePath * "/trials.jld2", "a+") do file
      ng = JLD2.Group(file, "trial_$(trial)")
      ng["finalTime"] = snapshot.time
      ng["front"] = snapshot.front
      ng["pIDs"] = getfield.(graph[sources], :ID_2)
      ng["mutID"] = getfield.(graph[sources], :ID_3)
    end

    # jldopen(save_path * "/data_sectors.jld2", "a+") do file
    #   ng = JLD2.Group(file, "trial_$trial")
    #   ng["sector_sizes"] = DataModels.lateralSectorSize(graph, cli.dims)
    # end

    jldopen(save_path * "/data_extras.jld2", "a+") do file
      ng = JLD2.Group(file, "trial_$trial")
      foreach(kv -> ng[String(first(kv))] = last(kv), unpack(dataModels))
    end
  end

  # .save snapshop of simulation
  if trial == 1
    if @isdefined(bndry)
      images(graph, width, height, cli, savePath, bndry)
    else
      images(graph, width, height, cli, savePath, missing)
    end
  end

  # .save ID_3 and update counts
  _resutls.heatmap .+= getfield.(graph, :ID_3)

  # .domain size scaling
  ξ = (0.0, 0.0)
  if cli.mutation > 0
    ξ, domain_hist = Observables.domainSizes(graph, cli.dims; nbins=10000)
    length(domain_hist) == length(_resutls.domainAreas) && (_resutls.domainAreas .+= domain_hist)
  end

  # .surviving numbers
  survFreq = Observables.survivalFreq(graph, sources, :ID_3)

  # .which hotspots contribute to the survival of mutants
  HotspotInfluences.contributing!(results, graph, cli)

  return round.((survFreq..., snapshot.time, snapshot.extinction, ξ...), digits=3)
end

# .delete item from namedtuple
delete(nt::NamedTuple{names}, keys::Vector{Symbol}) where {names} = NamedTuple{filter(x -> x ∉ keys, names)}(nt)
delete(nt::NamedTuple{names}, keys::Symbol) where {names} = NamedTuple{filter(x -> x != keys, names)}(nt)

# .command line options (input)
cli = Model.parseArgs()

# .set seed for reproducabilty
rngSeed = cli.rngSeed == 0 ? rand() : cli.rngSeed
seed!(rngSeed)

# .create containers for observales and associated functions
# dataModels = predefinedModels(animation = cli.animation)
dataModels = emptyModels(; animation=cli.animate)

# .create landscape of hotspots
env, htspts, cli = createEnvironment!(cli)

# .send command line parameters to stdout and mkdirs
cli = merge(cli, (rngSeed=rngSeed,))
savePath = setPath(cli)

# .dimensions
width = cli.width
height = cli.height

# .construct graphs
periodic = contains(cli.initial_type, "split") ? false : true
graph, cntns = hexGraph((width, height), periodic)

# .save all data in files
jldsave(savePath * "/inputs.jld2"; env=env, cli=cli, htspts=htspts, rngSeed=rngSeed)
jldsave(savePath * "/Opts.jld2"; delete(merge(cli, (rngSeed=rngSeed,)), :outputPath)...)

open(joinpath(savePath, "inputOpts.csv"), "w") do output
  opts = delete(merge(cli, (rngSeed=rngSeed,)), :outputPath)
  write(output, join(keys(opts), "\t") * "\n")
  write(output, join(values(opts), "\t") * "\n")
end

results = DataContainers.Results(
  zeros(Float32, size(graph)), zeros(Float32, 10000), zeros(UInt32, length(htspts)), htspts
)

open(joinpath(savePath, "table.csv"), "w") do output
  # .add header to table
  vars = "n_1,n_2,v_1,v_2,time,time_extinction,xi_m,xi_var"
  write(output, vars * "\n")

  for i in 1:(cli.numberTrials)
    data = main!(graph, i, dataModels, cli, env, results, savePath)
    write(output, join(data, ",") * "\n")
  end
end

#. save mutant spatial frequency data
if cli.heatmap
  open(Arrow.Writer, joinpath(savePath, "heatmap_ID3.arrow")) do file
    Arrow.write(file, (heatmap_ID3=results.heatmap,))
  end
end

# .no need to save this data if mutation rate is zero
if cli.mutation > 0
  open(Arrow.Writer, joinpath(savePath, "domain_sizes.arrow")) do file
    Arrow.write(file, (bin_counts=results.domainAreas,))
  end
end

# .save the most influential hotspots
open(Arrow.Writer, joinpath(savePath, "influentialHotspots.arrow")) do file
  Arrow.write(file, (influentialHotspots=results.influentialHotspots,))
end
