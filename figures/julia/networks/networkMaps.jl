using JLD2, FileIO, ArgParse, Base.Threads, Arrow
include("../..//src/base/model.jl")
include("../..//src/base/environment.jl")
include("../..//src/base/containers.jl")
include("../common/modeling.jl")

data_path = ""

#. interval sor sampling front
interval = 5
nfast = 1
NN = false

# .get experiment config
opts = namedtuple(load(joinpath(data_path, "Opts.jld2")))
opts = merge(opts, (ref_line=opts.height,))
(; width, height, radius, intensity, dims, ref_line) = opts

# .get environmental objects (hotspots)
if ispath(joinpath(data_path, "objects.jld2"))
  objs = load(joinpath(data_path, "objects.jld2"), "objs")
else
  objs, _env = load(joinpath(data_path, "inputs.jld2"), "htspts", "env")
end

#> NN hotspot graph
include("../calculations/OptimalPaths.jl")
htspts = filter(t -> last(t) <= opts.ref_line, unique!(objs))
points = transpose(hcat(Float64.(first.(htspts)), Float64.(last.(htspts))))

voronoiGraph = ContinuousOptPaths.voronoiTriangulation(points, opts)

# .nearest hotspot to each source (sink) node
sinks_to_htpts, sources_to_htspts = ContinuousOptPaths.sourceSinkConnections(
  voronoiGraph, opts; ref_line=opts.height, n_near=3, interval=interval
)

# .combine sources, sinks, and hotspots as nodes in the graph
g = ContinuousOptPaths.weightedGraph(voronoiGraph, sources_to_htspts, sinks_to_htpts, opts; NN=NN)

sourceIndices, sinkIndices = ContinuousOptPaths.edgeNodeIndices(
  sources_to_htspts, sinks_to_htpts, voronoiGraph.generators
)

gPositions = ContinuousOptPaths.graphPosisions(sources_to_htspts, sinks_to_htpts, voronoiGraph.generators, opts)

cutoff = 1.00001
optimalPathSets = @time ContinuousOptPaths.allOptimalPaths(
  sourceIndices, sinkIndices, opts, gPositions, g; n=200, cutoff=cutoff, floyd=false
)

# .reverse search, fastest paths from bottom to top instead of top to bottom
cutoff = 1.00001
optimalPathSets = @time ContinuousOptPaths.allOptimalPaths(
  sinkIndices, sourceIndices, opts, gPositions, g; n=200, cutoff=cutoff, floyd=false
)

#> snapshot of fastest paths and genetic lineages
include("../common/theme.jl")

shiftValues(x) = x == 1 ? 1 : 5

# .f_mt
hm = Arrow.Table(data_path * "/heatmap_ID3.arrow").heatmap_ID3
heatmap_ID3 = reshape(hm, opts.dims) ./ opts.numberTrials .- 1

# .influential heatmaps
influentialHotspots = convert(Vector{Float64}, Arrow.Table(joinpath(data_path, "influentialHotspots.arrow")).influentialHotspots)

influentialHotspots ./= opts.numberTrials

#. main plotting
begin
  lplot = false
  heatplot = true
  fig = Figure(;figure_padding=10)
  ax = Axis(fig[1, 1])
  hidedecorations!(ax)

  #. env and hotspots
  crange = (0.01, 0.4)
  heatplot && heatmap!(ax, heatmap_ID3; colormap=ColorSchemes.balance, colorrange=crange)
  Colorbar(fig[2, 1]; limits=crange, colormap=ColorSchemes.balance, vertical=false, flipaxis=false)
  rowgap!(fig.layout, 5)
  # heatmap!(ax, shiftValues.(reshape(env, opts.dims)); lowclip=:transparent, highclip=(:white, 0.65), colorrange=(2, 3))

  colors = get(cgrad(ColorSchemes.thermal), influentialHotspots, :extrema)
  scatter!(ax, getfield.(objs, 1), getfield.(objs, 2); color=colors, markersize=20)
  Colorbar(fig[1, 2]; limits=extrema(influentialHotspots), colormap=ColorSchemes.thermal)
  colgap!(fig.layout, 5)

  # .draw lineageMap by smoothing lineages with a 7x7 averaging kernel
  if lplot
    results::LineageTracing.Results = LineageTracing.lineageTraces(data_path, opts; full=false)
    smoothed_lineageMap = copy(results.lineageMap)
    for x in axes(smoothed_lineageMap, 1), y in axes(smoothed_lineageMap, 2)
      _xrange = clamp.((x - 3):(x + 3), 1, opts.width)
      _yrange = clamp.((y - 3):(y + 3), 1, opts.ref_line)
      smoothed_lineageMap[x, y] = sum(results.lineageMap[_xrange, _yrange]) / 49
    end

    smoothed_lineageMap ./= mean(smoothed_lineageMap[smoothed_lineageMap .> 0])
    heatmap!(ax, smoothed_lineageMap; lowclip=:transparent, colorrange=(0.03, 0.3), colormap=ColorSchemes.hot)
  end

  #. NN network paths
   for source_index in sinkIndices
    #for source_index in sourceIndices
    for optimalPath in optimalPathSets[source_index]
      xs = Float64[]
      ys = Float64[]
      for j in 2:2:length(optimalPath.pathPositions)
        j == firstindex(optimalPath.pathPositions) && continue
        (x1, y1) = optimalPath.pathPositions[j - 1]
        (x2, y2) = optimalPath.pathPositions[j]
        push!(xs, x1, x2)
        push!(ys, y1, y2)
      end
      lines!(ax, xs, ys; color=(:pink, 0.2), linewidth=1)
    end
  end

  ylims!(1, opts.ref_line)
  display(fig)
end