using CairoMakie, FileIO, LaTeXStrings

#doc conic section from Least Time Principle (at high intensity)
function conicSection(γ, cx, cy, R; resolution=100)
  θ = LinRange(0, π, resolution)
  x = @. cx + 2 * R * cos(θ) * sin(θ) / (γ - sin(θ))
  y = @. cy + 2 * R * sin(θ) * sin(θ) / (γ - sin(θ)) - R
  return x, y
end

function ellipseParameters(x, y)
  a = (maximum(x) - minimum(x)) / 2
  b = (maximum(y) - minimum(y)) / 2
  return a, b, hypot(a, b)
end

function centerEllipse(a, b, origin)
  x0 = origin[1] + a
  y0 = origin[2] + b
  return x0, y0
end

function fociEllipse(a, b, c, origin)
  x0, y0 = centerEllipse(a, b, origin)
  f1 = (x0 - c, y0)
  f2 = (x0 + c, y0)
  return f1, f2
end

#doc line segments approximating the tapering of bubbles
function edges(ν, x0, r, y0, cx, cy, resolution=100)
  y0 = y0 + cy - r
  x1 = LinRange(cx - x0, cx, resolution)
  x2 = LinRange(cx, cx + x0, resolution)

  l1 = @. (1 / (ν - 1)) * (x1 - (cx - 38)) + y0
  l2 = @. -(1 / (ν - 1)) * (x2 - (cx + 38)) + y0
  return l1, x1, l2, x2
end

#doc optics analogy of sector boundary angle (measured from y-axis)
function sectorBoundary(z, radius, intensity, selection, lx, ly, center)
#function sectorBoundary(z, radius, intensity, lx, ly, v_m)
  ν = 1 + intensity
  v_m = 1 - selection

  v_avg(z, r, ν) = v_m * z / (z - 2 * r + 2 * r / ν)
  θ = atan(2 * (v_avg(z, 1.15 * radius, ν) - 1) / (1 + v_avg(z, 1.15*radius, ν)))

  #r = sqrt(3) * radius
  r = radius

  #. center of the hotspot
  hx, hy = lx / 2, ly

  #xy_1 = @. hy * (sin(θ), cos(θ)) + (hx + 2 * r, 2 * r * sqrt(3))
  #xy_2 = @. hy * (-sin(θ), cos(θ)) + (hx - 2 * r, 2 * r * sqrt(3))
  #xy_10 = (hx + 2 * r, 2 * r * sqrt(3))
  #xy_20 = (hx - 2 * r, 2 * r * sqrt(3))

  xy_1 = @. hy * (sin(θ), cos(θ)) + (hx + 2 * r, 2 * r)
  xy_2 = @. hy * (-sin(θ), cos(θ)) + (hx - 2 * r, 2 * r)
  xy_10 = (hx + 2 * r, 2 * r)
  xy_20 = (hx - 2 * r, 2 * r)

  #.check if lines are intersecting
  m_1 = (xy_1[2] - xy_10[2]) / (xy_1[1] - xy_10[1])
  m_2 = (xy_2[2] - xy_20[2]) / (xy_2[1] - xy_20[1])
  if (m_1 == m_2)
    x_m = 0
  else
    x_m = (m_1 * xy_10[1] - m_2 * xy_20[1] + xy_20[2] - xy_10[2]) / (m_1 - m_2)
  end
  y_m = m_1 * (x_m - xy_10[1]) + xy_10[2]

  return xy_1, xy_10, xy_2, xy_20, (x_m, y_m)
end

function sortPaths(paths)
  sortby = x -> parse(Int, split(last(split(x, ",")), "_")[end])
  return sort(paths, by=sortby, rev=true)
end

dataAPI(path) =
  if isfile(joinpath(path, "inputs.jld2"))
    return load(joinpath(path, "inputs.jld2"), "cli", "env", "htspts")
  else
  end

#Z_c(r, ν, v) = (2r * ν - 2r) / (ν * (1 - v))
Z_c(r, ν, s) = (2 * r * ν) / (1 + ν) / s
shiftValues(x) = x == 1 ? 1 : 5

paths = readdir("/home/jgonzaleznunez/Projects/mutationWithLandscapes/workspace/simulations/transition_seq/"; join=true)
paths = filter(x -> isdir(x), paths)
paths = sortPaths(paths)

#let basePath = paths[15]
for basePath in paths
  mkpath(joinpath(basePath, "images"))

  #. extract data from file
  opts, env, htspts = dataAPI(basePath)

  heatmap_ID3 = reshape(load(basePath * "/heatmap_ID3.jld2", "heatmap_ID3"), opts.dims) ./ opts.numberTrials .- 1

  #. set up plot and axes with main data
  begin
    f_s = (550, 625)
    fig = Figure(; backgroundcolor=:white, size=f_s, figure_padding=20)
    ax = Axis(fig[2, 1])
    hidedecorations!(ax)

    if opts.separation < opts.height
      ax_2 = Axis(fig[2, 2]; yticklabelsvisible=true, xaxisposition=:bottom, yaxisposition=:right)
      ylims!(ax_2, (1, opts.height))
    end
    # xticks=round.([0, maximum(ordinate)], digits=3),
    # ticks=([0,0.5,1], ["0", "0.5", "1"]),

    kwargs = (ticksize=10, tickalign=1, ticklabelsize=15, colormap=:jet, vertical=false, flipaxis=true)
    Colorbar(fig[1, 1]; colorrange=(0, 0.75), kwargs...)
    # ticks=([0,0.5,1], ["0", "0.5", "1"]),
    # colorrange=extrema(heatmap_ID3),

    heatmap!(ax, heatmap_ID3; colormap=:jet, colorrange=(0, 0.75))
    heatmap!(
      ax, shiftValues.(reshape(env, opts.dims)); lowclip=:transparent, highclip=(:white, 0.95), colorrange=(2, 3)
    )

    rowgap!(fig.layout, 10)
    colgap!(fig.layout, 5)
    colsize!(fig.layout, 1, Relative(3 / 5))

    ylims!(ax, (1, opts.height))
    xlims!(ax, (1, opts.width))
  end

  if hasproperty(opts, :separation) && ~isnothing(opts.separation)
    # .conicSection of boundary  
    ν = 2 / (2 - opts.selection)
    cx, cy = first(htspts)
    xr, yr = conicSection(ν, cx, cy, opts.radius / 1.15; resolution=10000)
    scatter!(ax, xr, yr; color=:white, marker=:circle, markersize=5)

    # .stitch together the linear edges
    #xmax = maximum(xr .- cx)
    #left, x1, right, x2 = edges(ν, xmax, opts.radius, 100, cx, cy)
    #scatter!(ax, x1, left; color=:pink, marker=:dot, markersize=5)
    #scatter!(ax, x2, right; color=:pink, marker=:dot, markersize=5)
    #save(joinpath(basePath, "images", "sectorShapes.png"), fig)
  end

  if hasproperty(opts, :separation) && ~isnothing(opts.separation)
    #. if s > 0, then plot the critical separation for elliptical domains
    if opts.selection > 0
      ordinate = heatmap_ID3[cld(opts.width, 2), :]
      coordinate = collect(1:(opts.height))
      separation = 2 * sqrt(3) * opts.radius + opts.separation
      topElipse = maximum(yr)

      critS = Z_c(opts.radius, opts.intensity, opts.selection) + 2 * sqrt(3) * opts.radius
      hlines!(ax_2, critS; color=:green, label="Critical Separation")
    else
      ordinate = sum(heatmap_ID3; dims=1)[1, :] ./ size(heatmap_ID3, 1)
      coordinate = collect(1:(opts.height))
      separation = 0
      topElipse = 0
    end

    #. add separation lines to heatmap slice
    lines!(ax_2, ordinate, coordinate; color=:black)
    hlines!(ax_2, topElipse; color=:red, label="Ellipse Peak")
    hlines!(ax_2, separation; color=:blue, label="Separation", linewidth=5)

    save(joinpath(basePath, "images", "sectorShapes_$(opts.separation)_legendOff_boundariesOff.png"), fig)

    #. sector boundary angles
    p1, p10, p2, p20, xy_m = sectorBoundary(
      opts.separation,  opts.radius, opts.intensity, opts.selection, opts.width, opts.height, (0,0)
    )
    if xy_m[1] > 0 && xy_m[2] > 0
      p1 = xy_m
      p2 = xy_m
    end

    lines!(ax, [p10[1], p1[1]], [p10[2], p1[2]]; color=:green, linewidth=5)
    lines!(ax, [p20[1], p2[1]], [p20[2], p2[2]]; color=:green, linewidth=5)
  end

  save(joinpath(basePath, "images", "sectorShapes_$(opts.separation)_legendOff.png"), fig)
  axislegend(ax_2; merge=true, unique=true, position=:rt, fontsize=15)
  display(fig)
  #break
  save(joinpath(basePath, "images", "sectorShapes_$(opts.separation).png"), fig)
end