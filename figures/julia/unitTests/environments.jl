include("../../src/base/environment.jl")
using CairoMakie

begin
  opts = (radius=10, width=500, height=500, gap=10, separation=50, density=1, env_type="hex_grid", landscape=nothing)
  env, obs, opts = createEnvironment!(opts)

  fig = Figure()
  ax = Axis(fig[1, 1])

  heatmap!(
    ax,
    reshape(env, opts.width, opts.height);
    colormap=:grays,
    colorrange=(1.5, 1.6),
    lowclip=:transparent,
    highclip=:black
  )
  display(fig)
end