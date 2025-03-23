using FileIO, JLD2, StatsBase
include("../src/calculations/lineages.jl")
include("../src/calculations/OptimalPaths.jl")
include("common/theme.jl")
include("common/modeling.jl")
using .ContinuousOptPaths, DelaunayTriangulation

input = "/home/jgonzaleznunez/Projects/disorderedLandscapes/workspace/experiments/2024_18_02/XY:2000,1100_D:0.05_I:xxx_R:10_nEnvs:12_gap:0_intervals:0.0,2.0,8.0_trials:100_da/env_0/env_type_uniform,width_2000,height_1100,density_0.05,intensity_0.0,radius_10,numberTrials_100,numberSamples_50,rng_seed_1,rl_1000,da"

input = "/home/jgonzaleznunez/Projects/disorderedLandscapes/workspace/experiments/2024_18_02/XY:2000,1100_D:0.05_I:xxx_R:10_nEnvs:12_gap:0_intervals:0.0,2.0,8.0_trials:100_da/env_0/env_type_uniform,width_2000,height_1100,density_0.05,intensity_8.0,radius_10,numberTrials_100,numberSamples_50,rng_seed_1,rl_1000,da"

input = "/home/jgonzaleznunez/Projects/disorderedLandscapes/workspace/experiments/2024_25_02/XY:2000,1100_D:0.1_I:xxx_R:10_nEnvs:15_gap:0_intervals:0.0,2.0,8.0_trials:100_da/env_0/env_type_uniform,width_2000,height_1100,density_0.1,intensity_8.0,radius_10,numberTrials_100,numberSamples_50,rng_seed_1,rl_1000,da"
 
opts = namedtuple(load(joinpath(input, "Opts.jld2")))
htspts = load(joinpath(input, "objects.jld2"), "objs")

# >voronoi triangulation
htspts = filter(t -> last(t) <= opts.ref_line, unique!(htspts))
points = transpose(hcat(Float64.(first.(htspts)), Float64.(last.(htspts))))

voronoiGraph = ContinuousOptPaths.voronoiTriangulation(points, opts)

# >visualization of voronoi tessellation with nearest neighbors
function voronoiGraphImage(voronoiGraph, opts)
  fig = Figure(; size=(800, 800), figure_padding=30)
  ax = Axis(fig[1, 1]; title="Voronoi tessellation with graph", titlealign=:left)

  areas = get_area.(Ref(voronoiGraph.voronoi), 1:num_polygons(voronoiGraph.voronoi))
  colors = get(cgrad(ColorSchemes.grayC), areas, :extrema)

  voronoiplot!(ax, voronoiGraph.voronoi; show_generators=true, markercolor=:green, color=colors)

  for edge_candidate in voronoiGraph.edgeGenerators
    scatter!(ax, get_generator(voronoiGraph.voronoi, edge_candidate)...; color=:lightgreen, markersize=15)
  end

  for (gen_index, nbors) in voronoiGraph.connections
    for nbor in voronoiGraph.connections[gen_index]
      x1, y1 = get_generator(voronoiGraph.voronoi, gen_index)
      x2, y2 = get_generator(voronoiGraph.voronoi, nbor)
      abs(x2 - x1) > fld(opts.width, 2) && continue
      # task: color edges by distance, normalized by max distance between pairs
      lines!(ax, [x1, x2], [y1, y2]; color=:dodgerblue)
    end
  end

  limits!(ax, 0, opts.width, 0, opts.ref_line)
  display(fig)
end

voronoiGraphImage(voronoiGraph, opts)

# >nearest hotspot to each source (sink) node
sinks_to_htpts, sources_to_htspts = sourceSinkConnections(voronoiGraph, opts; ref_line=opts.ref_line, n_near=10)

# > combine sources, sinks, and hotspots as nodes in the graph
g = ContinuousOptPaths.weightedGraphComplete(voronoiGraph, sources_to_htspts, sinks_to_htpts, opts)

sourceIndices, sinkIndices = ContinuousOptPaths.edgeNodeIndices(voronoiGraph.generators, opts)
gPositions = ContinuousOptPaths.graphPosisions(voronoiGraph.generators, opts)
optimalPathSets = ContinuousOptPaths.allOptimalPaths(sourceIndices, sinkIndices, opts, gPositions, g; n=45, cutoff=1.15)

# > snapshot of optimal paths and lineage traces
function optimalPathsImage(sourceIndices, optimalPathSets, voronoiGraph, opts; vplot=false, lplot=false, pplot=false)
  fig = Figure(; size=(800, 600), figure_padding=20)
  ax = Axis(fig[1, 1]; backgroundcolor=:white)

  # .draw voronoi tessellation
  areas = get_area.(Ref(voronoiGraph.voronoi), 1:num_polygons(voronoiGraph.voronoi))
  colors = get(cgrad(ColorSchemes.grayC), areas, :extrema)

  vplot && voronoiplot!(ax, voronoiGraph.voronoi; show_generators=true, markercolor=:purple, color=colors)

  if lplot
    # .draw lineageMap by smoothing lineages with a 7x7 averaging kernel
    results::LineageTracing.Results = LineageTracing.lineageTraces(input, opts; full=false)
    smoothed_lineageMap = copy(results.lineageMap)
    for x in axes(smoothed_lineageMap, 1), y in axes(smoothed_lineageMap, 2)
      _xrange = clamp.((x - 3):(x + 3), 1, opts.width)
      _yrange = clamp.((y - 3):(y + 3), 1, opts.ref_line)
      smoothed_lineageMap[x, y] = sum(results.lineageMap[_xrange, _yrange]) / 49
    end

    smoothed_lineageMap ./= mean(smoothed_lineageMap[smoothed_lineageMap .> 0])
    heatmap!(ax, smoothed_lineageMap; lowclip=:transparent, colorrange=(0.01, 0.5), colormap=ColorSchemes.algae)
  end

  # .draw line segments between hotspots
  distances = Float64[]
  for source_index in sourceIndices
    push!(distances, getfield.(optimalPathSets[source_index], :pathDistance)...)
  end

  cmap = alphaColor(ColorSchemes.berlin, 0.8)
  colors = get(cgrad(cmap), distances, :extrema)

  if pplot
    for source_index in sourceIndices
      for optimalPath in optimalPathSets[source_index]
        # color by time deviation from minimum
        color = popfirst!(colors)

        for j in 2:2:length(optimalPath.pathPositions)
          (x1, y1) = optimalPath.pathPositions[j - 1]
          (x2, y2) = optimalPath.pathPositions[j]

          lines!(ax, [x1, x2], [y1, y2]; color=color, linewidth=5)
          # lines!(ax, [x1, x2], [y1, y2]; color=:blue, linewidth=5)
        end
      end
    end
  end

  hlines!(ax, 1000 - 100; color=:red)
  # for _visited_generatorIndex in _mask
  #   scatter!(ax, generators[_visited_generatorIndex]...; markersize=30, color=:pink)
  # end

  # .draw segments from top edges to hotspots
  # for i in axes(sources_to_htspts, 1)
  #   xposition, yposition = (Int64(i), opts.ref_line)
  #   for j in axes(sources_to_htspts, 2)
  #     centerX, centerY = edgeGencenters[sources_to_htspts[i, j]]
  #     lines!(ax, Int[xposition, centerX], Int[yposition, centerY]; color=(:blue, 0.2))
  #   end
  # end

  # # .draw segments from bottom edges to hotspots
  # for i in axes(sources_to_htspts, 1)
  #   xposition, yposition = (Int64(i), opts.ref_line)
  #   for j in axes(sources_to_htspts, 2)
  #     centerX, centerY = edgeGencenters[sources_to_htspts[i, j]]
  #     lines!(ax, Int[xposition, centerX], Int[yposition, centerY]; color=(:red, 0.2))
  #   end
  # end

  if pplot
    kwargs = (ticksize=15, tickalign=1, ticklabelsize=25, colormap=ColorSchemes.berlin, labelsize=25)
    # ticks = ([0, 0.5, 1], ["0", "0.5", "1"])
    colorrange = extrema(distances)
    Colorbar(fig[1, 2]; label="Path Travel Time", kwargs..., colorrange=colorrange)
  end

  if lplot
    kwargs = (
      ticksize=15,
      tickalign=1,
      ticklabelsize=25,
      colormap=ColorSchemes.algae,
      vertical=false,
      flipaxis=false,
      labelsize=25
    )
    # ticks = ([0, 0.5, 1], ["0", "0.5", "1"])
    colorrange = (0.01, 0.5)
    Colorbar(fig[2, 1]; label="Lineage Visitation Frequency", kwargs..., colorrange=colorrange)
  end

  scalebar!(ax, (110, 100); width=200, height=20, color=:red, fontsize=35)
  limits!(ax, 1, opts.width, 1, opts.ref_line)
  colgap!(fig.layout, 1)
  rowgap!(fig.layout, 1)
  hidedecorations!(ax)
  display(fig)
end

optimalPathsImage(sourceIndices, optimalPathSets, voronoiGraph, opts; lplot=false, pplot=true, vplot=true)

optimalPathsImage(sourceIndices, optimalPathSets, voronoiGraph, opts; lplot=true, pplot=true, vplot=false)

# >what is the MSD of the optimal paths
optimalPathResults = ContinuousOptPaths.optimalPathMSD(optimalPathSets, opts.ref_line, opts.width)

fig = Figure(; size=(800, 600), figure_padding=20)
ax = Axis(fig[1, 1]; backgroundcolor=:white)

areas = get_area.(Ref(voronoiGraph.voronoi), 1:num_polygons(voronoiGraph.voronoi))
colors = get(cgrad(ColorSchemes.grayC), areas, :extrema)
voronoiplot!(ax, voronoiGraph.voronoi; show_generators=true, markercolor=:purple, color=colors)

begin
  source_index = rand(sourceIndices)
  for optimalPath in optimalPathSets[source_index]
    for j in 2:2:length(optimalPath.pathPositions)
      (x1, y1) = optimalPath.pathPositions[j - 1]
      (x2, y2) = optimalPath.pathPositions[j]
      lines!(ax, [x1, x2], [y1, y2]; color=:red, linewidth=5)
    end
    println(first(optimalPath.pathPositions)[1] - last(optimalPath.pathPositions)[1])
    println(first(optimalPath.pathPositions), " ", last(optimalPath.pathPositions))
    break
  end
  display(fig)
end

_msd, _msd_2 = let
  _msd = zeros(Float64, 1000)
  _counts = zeros(UInt32, 1000)

  _msd_2 = zeros(Float64, 1000)
  _counts_2 = zeros(UInt32, 1000)

  # task: set x reference point
  
end

x = collect(1:1000)
# fig, ax = scatter(x[_msd .> 0], _msd[_msd .> 0], color=:white)
fig, ax = scatter(x[2:end], _msd_2[2:end], color=:black)
# scatter!(ax, x, _msd_2, color=:black)
# scatter!(ax, x[_msd .> 0], _msd[_msd .> 0], color=:red)
# scatter!(ax, x[_msd .> 0], _msd_2[_msd .> 0], color=:black)
ax.xscale = log10
ax.yscale = log10
display(fig)

x = collect(1:(opts.ref_line))
y = optimalPathResults.MSD

fit = FitModels.linearfit(x, y; xlbound=10, xhbound=600, scale=log10)

maskedX, maskedY = FilterTools.mask(x, y; xlbound=10, xhbound=600, scale=log10)

x, y = FilterTools.mask(x, y; xlbound=1, scale=log10)

let
  fig, ax = scatter(x, y)
  scatter!(ax, maskedX, maskedY; color=:red)

  # fig, ax = scatter((getfield.(unique(optPath.pathPositions), 1) .- 682)[2:end])
  ax.xlabel = "vertical distance from front"
  ax.ylabel = "Lateral MSD"
  ax.xscale = log10
  ax.yscale = log10
  # limits!(ax, (1, 1200), (1, 300))
  display(fig)
end

# >visited hotspots by optimal paths
visitedObjectsIndices = ContinuousOptPaths.objectVisits(voronoiGraph, sourceIndices, optimalPathSets)

# >is there a correlation between cell area and visitation
reached = zeros(Int, length(kdtree.data))
cell_areas = get_area.(Ref(_voronoi), collect(keys(generators)))
_mask = _generatorIndices[findall(>(0.05), reached ./ maximum(reached))]

# >how does path length compare to shortest line distance
for source_index in sourceIndices
  push!(distances, getfield.(optimalPathSets[source_index], :pathDistance)...)
end

# >what is the average neighbor cell sizes (or distance) compared to visited cell