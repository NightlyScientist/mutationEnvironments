include("../../src/base/environment.jl")
include("../../src/base/model.jl")
include("../../src/base/trees.jl")
include("../../src/base/dataModels.jl")
include("../common/theme.jl")
include("../common/modeling.jl")

# .global parameters
dimensions = 500, 500
regionBounds = 1, dimensions[2]
density = 0.05
intensity = 3
radius = 50
stopline = 17
env_type = "centered_circle"
gap = 5
selection = 0.03
mutation = 0.0001
standing_variation = false

# .base path for figures and videos
basePath = "/home/jgonzaleznunez/Projects/mutationWithLandscapes/workspace/media/single_hotspot/"
mkpath(basePath)

cli = (
  height=dimensions[2],
  width=dimensions[1],
  env_type=env_type,
  radius=radius,
  density=density,
  gap=gap,
  dims=dimensions,
  numberSamples=50,
  selection=selection,
  intensity=intensity,
  mutation=mutation,
  standing_variation=standing_variation,
  singleMutant=true,
  landscape=nothing,
  initial_type="alt",
)

customTheme!(20)
label = :ID_3
cmap = ColorScheme([RGB(1, 0, 0), RGB(1, 1, 0)])

begin
  # .generate enviroment of hotspots
  env, htspts, cli = createEnvironment!(cli)

  # .do simulation
  dataModels = emptyModels(; animation=true)
  graph, cntns = Model.hexGraph(dimensions, true)
  snapshot = Model.simulate!(graph, cntns, cli, dataModels, env)
  sources = snapshot.front

  # .determine the phylo of the periphery population 
  phylo, _ = GenealogicalTree.genealogy(graph, sources)

  fig = Figure(; size=(500, 500 * sqrt(3) / 2), backgroundcolor=:white)
  ax = Axis(fig[1, 1])
  hidedecorations!(ax)

  # .grab snapshots from animation sequence and take the last one
  snapshots = dataModels[:animation].data

  i = lastindex(snapshots)
  time, ids_2, ids_3, lineages = snapshots[i]

  labelField = reshape(getfield.(graph, label), cli.dims)

  if label == :ID_1
    # .don't draw hotspots if intensity is zero
    remapping = cli.intensity == 0.0 ? Int64[1, 1, 2, 2] : Int64[1, 3, 2, 3]
    labelField = map(t -> remapping[t], labelField)
    colorrange = (1, 2)
  end

  heatmap!(ax, labelField; colormap=cmap, colorrange=(1, 2), lowclip=:white, highclip=:transparent)

  # .hotspots
  heatmap!(ax, reshape(env, cli.dims); colorrange=(1.1, 1.2), lowclip=:transparent, highclip=(:black, 0.7))

  # .save images to path
  # save(joinpath(basePath, "single_hotspot_$i.png"), fig)

  display(fig)
end