using JLD2, FileIO, Arrow
import StatsBase.Random: seed!
include("../../calculations/scaling.jl")
include("../../base/dataModels.jl")
include("../../alt_models/oneDimensional/model.jl")
include("../../base/environment.jl")
include("../../base/trees.jl")

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

function images(graph, env, width, height, cli, savePath, bndry=missing)
  for label in [:ID_1, :ID_2, :ID_3, :ID_4]
    labelField = reshape(getfield.(graph, label), (width, height))

    fig = Figure(; size=(500, 0.87 * 500 * height / width))
    ax = Axis(fig[1, 1])

    if (label == :ID_2 || label == :ID_4)
      colorrange = (1, width)
      cmap = tab10
    elseif label == :ID_1
      colorrange = (1, 6)
      cmap = ColorSchemes.roma
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
    hidedecorations!(ax)
    save(joinpath(savePath, "snapshots_$(label).png"), fig)
  end
end

function main!(graph, trial, dataModels, cli, env, results, data_path)
  periodic = contains(cli.initial_type, "split") ? false : true

  # .reset models
  resetModels!(dataModels)

  # .reset graph properties
  resetGraph!(graph)

  snapshot = simulate!(graph, cli, dataModels, env; periodic=periodic)
  sources = Model.graphSources(cli)

  if cli.detailed_analytics
    # .determine the phylo of the periphery population 
    phylo, branchPoints, branchTimes = GenealogicalTree.genealogy(graph, sources)
    geneticLabels = getfield.(graph[sources], :ID_2)

    jldopen(data_path * "/data_phylo.jld2", "a+") do file
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

    jldopen(data_path * "/data_sectors.jld2", "a+") do file
      ng = JLD2.Group(file, "trial_$trial")
      ng["sector_sizes"] = DataModels.lateralSectorSize(graph, opts.dims)
    end

    jldopen(data_path * "/data_extras.jld2", "a+") do file
      ng = JLD2.Group(file, "trial_$trial")
      foreach(kv -> ng[String(first(kv))] = last(kv), DataModels.unpack(dataModels))
    end
  end

  # .trace boundary 
  if contains(cli.initial_type, "split") && cli.standing_variation
    bndry = Observables.traceBoundary(graph, width, height; window=20, delim=3, stop_index=snapshot.stop_index)
    if trial == 1
      results.sectorMV = zeros(Float32, length(bndry), 3)
      results.sectorMV[:, 1] .= 1
      results.sectorMV[:, 2] .= bndry
    else
      Observables.onlineUpdates!(results.sectorMV, bndry)
    end
  end

  #. track mutation density with expansion distance
  results.mutationalFreq .+= Observables.mutantFrequency(graph, width, height)

  # .save snapshop of simulation
  if trial == 1
    if @isdefined(bndry)
      images(graph, env, width, height, cli, savePath, bndry)
    else
      images(graph, env, width, height, cli, savePath, missing)
    end
  end

  # .save ID_3 and update counts
  results.heatmap .+= getfield.(graph, :ID_3)

  # .domain size scaling
  ξ = (0.0, 0.0)
  if cli.mutation > 0
    ξ, domain_hist = Observables.domainSizes(graph, cli.dims; nbins=10000)
    length(domain_hist) == length(results.domainAreas) && (results.domainAreas .+= domain_hist)
  end

  # .surviving numbers
  survFreq = Observables.survivalFreq(graph, sources, :ID_3; species=3)

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
dataModels = emptyModels(; animation=cli.animate)

# .create landscape of hotspots
env, htspts, cli = createEnvironment!(cli)

# .assert proper bounds of input parameters, c.g. Castillo 2020
if ~(0 <= cli.selection < 1 && 0 <= cli.compensation <= cli.selection)
  println("check input parameters\n")
  @warn "check input parameters\n"
  exit()
end

# .send command line parameters to stdout and mkdirs
cli = merge(cli, (rngSeed=rngSeed,))
savePath = setPath(cli)

# .dimensions
width = cli.width
height = cli.height

# .construct graphs
graph = hexGraph((width, height))

# .save all data in files
jldsave(savePath * "/inputs.jld2"; env=env, cli=cli, htspts=htspts, rngSeed=rngSeed)
jldsave(savePath * "/objects.jld2"; objs=htspts)
jldsave(savePath * "/Opts.jld2"; delete(merge(cli, (rngSeed=rngSeed,)), :outputPath)...)

open(joinpath(savePath, "inputOpts.csv"), "w") do output
  opts = delete(merge(cli, (rngSeed=rngSeed,)), :outputPath)
  write(output, join(keys(opts), "\t") * "\n")
  write(output, join(values(opts), "\t") * "\n")
end

#. data used across many independent runs
results = DataContainers.Results()
results.heatmap = zeros(Float32, size(graph))
results.domainAreas = zeros(Float32, 10000)
results.htspts = htspts
results.mutationalFreq = zeros(Float32, height)

open(joinpath(savePath, "table.csv"), "w") do output
  # .add header to table
  vars = "n_1,n_2,n_3,v_1,v_2,v_3,time,time_extinction,xi_m,xi_var"
  write(output, vars * "\n")

  for i in 1:(cli.numberTrials)
    data = main!(graph, i, dataModels, cli, env, results, savePath)
    write(output, join(data, ",") * "\n")
  end
end

if cli.heatmap
  open(Arrow.Writer, joinpath(savePath, "heatmap_ID3.arrow")) do file
    Arrow.write(file, (heatmap_ID3=results.heatmap,))
  end
end

#. save sector boundary data
ismissing(results.sectorMV) || open(Arrow.Writer, joinpath(savePath, "boundary.arrow")) do file
  _sectorPos, _, _sectorVar = Observables.retrieveOnlineUpdates(results.sectorMV)
  Arrow.write(file, (sectorPos=_sectorPos, sectorVar=_sectorVar))
end

#. save mutation frequency
open(Arrow.Writer, joinpath(savePath, "mutationalFreq.arrow")) do file
  _mf = results.mutationalFreq ./ cli.numberTrials
  Arrow.write(file, (mutationalFreq=_mf,))
end

# .save mutation bubble (sector) data
if cli.mutation > 0
  open(Arrow.Writer, joinpath(savePath, "domain_sizes.arrow")) do file
    Arrow.write(file, (bin_counts=results.domainAreas,))
  end
end