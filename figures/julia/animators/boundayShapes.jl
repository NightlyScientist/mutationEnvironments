using CairoMakie, FileIO, LaTeXStrings

function ellipse(ν, cx, cy, R, resolution=100)
  θ = LinRange(0, π, resolution)
  x = @. cx + 2 * R * sin(θ) * cos(θ) / (ν - sin(θ))
  y = @. cy .+ 2 * R .* sin.(θ) .* sin.(θ) / (ν - sin(θ)) - R
  return x, y
end

function edges(ν, x0, r, y0, cx, cy, resolution=100)
  y0 = y0 + cy - r
  x1 = LinRange(cx - x0, cx, resolution)
  x2 = LinRange(cx, cx + x0, resolution)

  left = @. (1 / (ν - 1)) * (x1 - (cx - 38)) + y0
  right = @. (-1 / (ν - 1)) * (x2 - (cx + 38)) + y0
  return left, x1, right, x2
end

function slopes(s, radius, n, lx, ly, vm)
  # r = radius
  r = ceil(sqrt(3) * radius)
  v_avg(s, r, n) = vm * s / (s - 2 * r + 2 * r / n)

  θ = atan(2 * (v_avg(s, radius, n) - 1) / (1 + v_avg(s, radius, n)))

  xy_1 = ly .* (sin(θ), cos(θ)) .+ (lx / 2 + 2 * r, 2 * r * sqrt(3))
  xy_10 = (lx / 2 + 2 * r, 2 * r * sqrt(3))
  xy_2 = ly .* (-sin(θ), cos(θ)) .+ (lx / 2 - 2 * r, 2 * r * sqrt(3))
  xy_20 = (lx / 2 - 2 * r, 2 * r * sqrt(3))

  return xy_1, xy_10, xy_2, xy_20
end

separation_crit(r, n, v) = (2r * n - 2r) / (n * (1 - v))

shiftValues(x) = x == 1 ? 1 : 5
# ****

inputs = readdir(
  "/home/jgonzaleznunez/Projects/mutationWithLandscapes/workspace/experiments/boundary_shape/XY:500,500_S:0.1_I:4.0_m:0.0_d:0.07_r:10_intervals:20.0,10.0,400.0_trials:500_EV:circle";
  join=true
)

img_dir = "/home/jgonzaleznunez/Projects/mutationWithLandscapes/workspace/animations"

sorted_inputs = let opt = :separation, inputs = inputs
  opt_values = []
  inputs = filter(t -> ~contains(t, "logs"), inputs)
  for input in inputs
    opts = load(input * "/inputs.jld2", "cli")
    push!(opt_values, get(opts, opt, 0))
  end
  println(opt_values[sortperm(opt_values; rev=true)])
  inputs[sortperm(opt_values; rev=true)]
end

begin
  fig = Figure(; backgroundcolor=:white, resolution=(600, 650))
  # fig = Figure(resolution=(600, 600 * 1.00))
  ax = Axis(
    fig[2, 1]; leftspinevisible=false, bottomspinevisible=false, yticklabelsvisible=false, xticklabelsvisible=false
  )

  ax_2 = Axis(fig[2, 2]; yticklabelsvisible=true, xaxisposition=:bottom, yaxisposition=:right)
  # xticks=round.([0, maximum(ordinate)], digits=3),
  # ticks=([0,0.5,1], ["0", "0.5", "1"]),

  kwargs = (ticksize=15, tickalign=1, ticklabelsize=25, colormap=:jet, vertical=false, flipaxis=true, labelsize=35)
  Colorbar(fig[1, 1]; label=L"Mutant Frequency $M(x,y)$", colorrange=(0, 0.7), kwargs...)
  # ticks=([0,0.5,1], ["0", "0.5", "1"]),
  # colorrange=extrema(heatmap_ID3),

  xlims!(ax_2, (0, 0.7))

  record(fig, joinpath(img_dir, "boundary_shapes.mp4"); framerate=2) do io
    for input in sorted_inputs
      opts, env, htspts = load(joinpath(input, "inputs.jld2"), "cli", "env", "htspts")

      heatmap_ID3 = reshape(load(input * "/heatmap_ID3.jld2", "heatmap_ID3"), opts.dims) ./ opts.numberTrials .- 1

      empty!(ax)
      empty!(ax_2)

      heatmap!(ax, heatmap_ID3; colormap=:jet)
      heatmap!(
        ax, shiftValues.(reshape(env, opts.dims)); lowclip=:transparent, highclip=(:white, 0.65), colorrange=(2, 3)
      )

      if hasproperty(opts, :separation) && ~isnothing(opts.separation)
        # .ellipse of boundary  
        ν = 1 / (1 - opts.selection)
        cx, cy = first(htspts)
        xr, yr = ellipse(ν, cx, cy, opts.radius)
        xmax = maximum(xr .- cx)
        scatter!(ax, xr, yr; color=:white, markerstyle=:dot, markersize=10)

        # .stitch together the linear edges
        left, x1, right, x2 = edges(ν, xmax, opts.radius, 100, cx, cy)
        scatter!(ax, x1, left; color=:pink, markerstyle=:dot, markersize=5)
        scatter!(ax, x2, right; color=:pink, markerstyle=:dot, markersize=5)
      end

      colgap!(fig.layout, 2)
      rowgap!(fig.layout, 2)

      # >mutFrequency
      if hasproperty(opts, :separation) && ~isnothing(opts.separation)
        ordinate = heatmap_ID3[cld(opts.width, 2), :]
        coordinate = collect(1:(opts.height))
        separation = 2 * sqrt(3) * opts.radius + opts.separation
        topElipse = maximum(yr)
      else
        ordinate = sum(heatmap_ID3; dims=1)[1, :] ./ size(heatmap_ID3, 1)
        coordinate = collect(1:(opts.height))
        separation = 0
        topElipse = 0
      end

      lines!(ax_2, ordinate, coordinate; color=:black)
      hlines!(ax_2, topElipse; color=:red, label="Ellipse Peak")
      hlines!(ax_2, separation; color=:blue, label="Separation", linewidth=5)

      critS = separation_crit(opts.radius, opts.intensity + 1, 1 - opts.selection)
      hlines!(ax_2, 2 * sqrt(3) * opts.radius + critS; color=:green, label="Critical Separation")

      ylims!(ax, (1, opts.height))
      xlims!(ax, (1, opts.height))
      ylims!(ax_2, (1, opts.height))
      # xlims!(ax_2, (0, maximum(ordinate)))

      p1, p10, p2, p20 = slopes(
        opts.separation, opts.radius, opts.intensity + 1, opts.width, opts.height, 1 - opts.selection
      )
      lines!(ax, [p10[1], p1[1]], [p10[2], p1[2]]; color=:green, linewidth=3)
      lines!(ax, [p20[1], p2[1]], [p20[2], p2[2]]; color=:green, linewidth=3)

      colgap!(fig.layout, 0)
      colsize!(fig.layout, 1, Relative(4 / 5))
      Legend(fig[1, 2], ax_2, ; nbanks=1)
      # display(fig)

      recordframe!(io)  # record a new frame
    end
  end
end
