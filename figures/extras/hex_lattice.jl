include("../../src/base/environment.jl")
include("../../src/base/model.jl")
include("../../src/base/trees.jl")
include("../../src/base/dataModels.jl")
include("../common/theme.jl")
include("../common/modeling.jl")

# .global parameters
dimensions = 20, 20
regionBounds = 1, dimensions[2]
density = 0.05
intensity = 3
radius = 2
stopline = 17
env_type = "uniform"
gap = 5
selection = 0.1
mutation = 0.1
standing_variation = false

# .base path for figures and videos
basePath = "/home/jgonzaleznunez/Projects/mutationWithLandscapes/workspace/media/hex_lattice/"
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
  singleMutant=false
)

env, htspts, cli = createEnvironment!(cli)

# .do simulation
dataModels = emptyModels(; animation=true)
graph, cntns = Model.hexGraph(dimensions, true)
snapshot = Model.simulate!(graph, cntns, cli, dataModels, env)
sources = snapshot.front

# .determine the phylo of the periphery population 
phylo, _ = GenealogicalTree.genealogy(graph, sources)

customTheme!(20)

x, y = reconstructCoordinates(collect(1:length(graph)), cli.width, cli.height)
gy = (y .- 1) .* 1.5 .+ 1.5
gx = (x .- 0.5 * isodd.(y)) .* (sqrt(3))

let canvasLX = 450, ms = canvasLX / 16
  fig = Figure(; size=(canvasLX, canvasLX * sqrt(3) / 2), backgroundcolor=:white)
  ax = Axis(fig[1, 1])
  hidedecorations!(ax)

  snapshots = dataModels[:animation].data

  # record(fig, joinpath(videoPath, "growth.mp4"); framerate=12) do io
  for i in eachindex(snapshots)
    time, ids_2, ids_3, lineages = snapshots[i]

    # .empty previous content from axis
    # empty!(ax)

    # .grid as hexagons (black outline)
    scatter!(ax, gx, gy; marker=:hexagon, markersize=ms + 0.5, color=:black)

    # .grid as hexagons (white fill)
    scatter!(ax, gx, gy; marker=:hexagon, markersize=ms - 0.5, color=(:white, 1.0))

    # .color initial population
    color = collect(1:dimensions[1])
    # .cylic colormap
    cmap = ColorScheme([RGB(1, 1, 1), RGB(1, 0, 0), RGB(0, 0, 0)])
    alphaColors = alphaColor(cmap, 1; ncolors=cli.dims[1])
    scatter!(
      ax,
      gx[color],
      gy[color];
      marker=:hexagon,
      markersize=ms - 10,
      color=getfield.(graph[color], :ID_3),
      colormap=alphaColors
    )

    # .color population
    # times = getfield.(graph, :time)
    # fltr = times .<= maximum(times[snapshot.front])
    # fltr[color] .= false
    # color = getfield.(graph, :ID_3)
    cmap = [RGB(1, 1, 1), RGB(1, 0, 0), RGB(1, 1, 0), RGB(0, 1, 0)]
    # colors = map(t -> cmap[t], color)
    # scatter!(ax, gx[fltr], gy[fltr]; marker=:hexagon, markersize=ms - 1.5, color=colors[fltr])
    colors = map(t -> cmap[t], ids_3 .+ 1)
    scatter!(ax, gx[21:end], gy[21:end]; marker=:hexagon, markersize=ms - 1.5, color=colors[21:end])

    # .highlight the front
    # scatter!(ax, gx[snapshot.front], gy[snapshot.front]; marker=:hexagon, markersize=ms, color=(:blue, 0.7))

    # .hotspots
    scatter!(ax, gx[env .== 2], gy[env .== 2]; marker=:hexagon, markersize=ms - 1, color=(:black, 0.4))

    # .draw lineages as sequence of segments
    # for source in sources
    #   current = source
    #   idxs = Int64[]

    #   while current != 0
    #     push!(idxs, current)
    #     current = phylo[current]
    #   end
    #   lines!(ax, gx[idxs], gy[idxs]; color=(:green, 1), linewidth=5)
    # end

    save(joinpath(basePath, "hexlattice_$i.png"), fig)
    # recordframe!(io)  # record a new frame
  end
  # end
  display(fig)
end