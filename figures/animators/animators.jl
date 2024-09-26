include("../../src/base/model.jl")
include("../common/modeling.jl")
include("../common/theme.jl")

using JLD2, CairoMakie, ColorSchemes, Colors
using ArgParse
import .Model: gNodes, buildMap
import ColorSchemes: tab10, ColorScheme, berlin25
import Colors: RGBA

addArgs!(sts, name; kwargs...) = add_arg_table!(sts, "--$(name)", Dict(kwargs))

function ParseArgs()
  sts = ArgParseSettings()
  addArgs!(sts, "input"; required=true)
  addArgs!(sts, "output"; required=true)
  addArgs!(sts, "overwrite"; action=:store_true)
  addArgs!(sts, "fps"; default=14, arg_type=Int)
  addArgs!(sts, "lineages"; action=:store_true)
  addArgs!(sts, "hotspots"; action=:store_true)
  addArgs!(sts, "colors"; action=:store_true)
  addArgs!(sts, "mutations"; action=:store_true)
  return namedtuple(parse_args(sts))
end

cli = ParseArgs()
@info "cli:" cli

if ~cli.lineages && ~cli.colors && ~cli.mutations
  @warn "need at least one of: colors, lineages"
  exit()
end

@info "saving recording to:" cli.output

fn = let
  s = splitdir(cli.input)
  last(s) == "" ? last(splitdir(first(s))) : last(s)
end

@info "using data from: " fn
mkpath(cli.output * "/" * fn)

let cmdline = cli
  env, cli = load(joinpath(cmdline.input, "inputs.jld2"), "env", "cli")

  file = jldopen(joinpath(cmdline.input, "data_extras.jld2"), "r")

  cli.animate == false && return println("no animation found")

  # theme and color
  set_theme!(; figure_padding=00)
  colors = alphaColor(ColorSchemes.grays, 0.3; ncolors=3)
  greens = alphaColor(ColorSchemes.Greens, 0.5; ncolors=2)

  cmap = ColorScheme([RGB(1, 0, 0), RGB(1, 1, 0)])

  # create figure
  fig = Figure(; size=cli.dims .* (√(3), 1.5), backgroundcolor=:transparent)
  ax = Axis(fig[1, 1]; backgroundcolor=:white)

  # set image limits
  hidedecorations!(ax)
  limits!(ax, 1, cli.dims[1], 1, cli.dims[2])

  # used for updating lineages
  tmp = zeros(Int64, size(env))

  shiftValues(x) = x == 0 ? 0 : 5

  record(fig, joinpath(cmdline.output, fn, "growth.mp4"); framerate=cmdline.fps) do io
    for trial in FilterTools.mask(keys(file), "trial")
      # create object for snapshots: list->tuple(growth, lineages)
      snapshots = file[trial]["animation"]

      for i in eachindex(snapshots)
        # empty previous content from axis
        empty!(ax)
        time, ids_2, ids_3, lineages = snapshots[i]

        # hack: needed to comport with existing variables
        ids = ids_2

        # id colors
        if cmdline.colors
          heatmap!(ax, reshape(ids, cli.dims); colormap=:hsv, colorrange=(1, cli.dims[1]), lowclip=:transparent)
        elseif cmdline.mutations
          heatmap!(ax, reshape(ids_3, cli.dims); colormap=cmap, colorrange=(1, 2), lowclip=:transparent)
        else
          heatmap!(
            ax, shiftValues.(reshape(ids, cli.dims)); colorrange=(1, 2), lowclip=:transparent, highclip=(:green, 0.2)
          )
        end

        # lineages
        if cmdline.lineages && length(lineages) > 0
          tmp .= 0
          tmp[lineages] .= 3
          color = cmdline.colors ? :black : :black
          heatmap!(ax, reshape(tmp, cli.dims); colorrange=(1, 2), lowclip=:transparent, highclip=color)
        end

        # hotspots
        if cmdline.hotspots
          heatmap!(ax, reshape(env, cli.dims); colorrange=(2, 3), lowclip=:transparent, colormap=colors)
        end
        recordframe!(io)  # record a new frame
      end
    end
  end

  close(file)
end
