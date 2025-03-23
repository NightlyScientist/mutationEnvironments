using JLD2, FileIO
include("common/theme.jl")
include("common/binning.jl")
include("common/modeling.jl")
include("../src/calculations/lineages.jl")
include("../src/base/environment.jl")

function saveImage(fig, path, imgPath, name)
  imgFolderPath = last(splitdir(path))
  img_path = joinpath(imgPath, imgFolderPath)
  figPath = joinpath(img_path, name)
  mkpath(img_path)
  save(figPath, fig)
  return figPath
end

# doc: add heat map of surviving ancestors below lineage map
function addSurvivalFreqImg!(fig, result, cmap)
  sub_axis = Axis(fig[3, 1]; backgroundcolor=:white, xlabelsize=24, ylabelsize=24)

  freq_ancestors = result.ancestorCounts ./ sum(result.ancestorCounts)
  img = reshape(repeat(freq_ancestors, 2), (length(freq_ancestors), 2))
  heatmap!(sub_axis, img; lowclip=:white, colormap=cmap)

  hidedecorations!(sub_axis; grid=true)

  kwargs = (ticksize=5, ticklabelsize=20, vertical=false, flipaxis=false, labelsize=24, label="Survival Probability")
  Colorbar(fig[4, 1]; colormap=cmap, ticks=([0, 0.5, 1], ["0", "0.5", "1"]), kwargs...)
  xlims!(sub_axis, (1, length(freq_ancestors)))
  rowsize!(fig.layout, 3, Fixed(25))
  rowsize!(fig.layout, 4, Fixed(5))
  rowgap!(fig.layout, 1)
  return nothing
end

# doc: add main heat map of lineage positions
function addLineagemapImg!(fig, result, cmap, opts)
  ax = Axis(fig[2, 1]; backgroundcolor=:white)

  # task: normalize to something better
  hm = result.lineageMap ./ maximum(result.lineageMap)
  colorrange = (0.001, 0.5 * maximum(hm))
  h = heatmap!(ax, hm; lowclip=:white, colormap=cmap, colorrange=colorrange)

  kwargs = (ticksize=10, ticklabelsize=20, tickalign=1, labelsize=24, flipaxis=true, vertical=false)
  Colorbar(fig[1, 1], h; label="Lineage Visitation Frequency", kwargs...)

  colors = alphaColor(ColorSchemes.grays, 0.3; ncolors=3)
  if opts.intensity > 0
    objs = load(joinpath(data_path, "objects.jld2"), "objs")
    env = first(applyObstacles!(objs, opts.radius, opts.width, opts.height))
    env = reshape(env, opts.dims)
    heatmap!(ax, env; colorrange=(2, 3), lowclip=:transparent, colormap=colors)
  end

  rowgap!(fig.layout, 1)
  colgap!(fig.layout, 1)
  rowsize!(fig.layout, 2, Fixed(300))
  rowsize!(fig.layout, 1, Fixed(10))
  hidedecorations!(ax)
  scalebar!(ax, (150, 50); width=200, height=3, offset=50, color=:black)
  limits!(ax, 1, opts.width, 1, opts.ref_line)
  return nothing
end
# >collect data files
# task: convert to dataframes
data_path = "/home/jgonzaleznunez/Projects/disorderedLandscapes/workspace/experiments/2024_25_02/XY:2000,1100_D:0.1_I:xxx_R:10_nEnvs:15_gap:0_intervals:0.0,2.0,8.0_trials:100_da/env_0/env_type_uniform,width_2000,height_1100,density_0.1,intensity_8.0,radius_10,numberTrials_100,numberSamples_50,rng_seed_1,rl_1000,da"

_dataFileName = "data_phylo.jld2"
_fieldNames = ["end_time", "source_labels", "phylogeny", "branch_points", "branch_times"]

opts = namedtuple(FileIO.load(joinpath(data_path, "Opts.jld2")))
results::LineageTracing.Results = LineageTracing.lineageTraces(data_path, opts)

customTheme!(22)
cmap = alphaColor(ColorSchemes.dense, 1.0)
cmap_2 = alphaColor(ColorSchemes.tempo, 1.0)

# >pinning heatmap and ancestor Survival Frequencies
fig = Figure(; backgroundcolor=:white, size=(600, 600 * 0.87),)

# .lineage heat map of Visitation Frequency
addLineagemapImg!(fig, results, cmap, opts)

# .surviving ancestors heatmap
addSurvivalFreqImg!(fig, results, cmap_2)
display(fig)