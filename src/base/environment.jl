using FileIO, JLD2, DelimitedFiles

function createEnvironment!(opts)
  upnt(c, k, v) = merge(c, (Symbol(k) => v,))

  opts = upnt(opts, "dims", (opts.width, opts.height))

  if ~isnothing(opts.landscape)
    ispath(opts.landscape) || (@info "could not find environment.txt in path."; exit())
    objs = readdlm(opts.landscape, ',', Int, '\n')
    obs = tuple.(eachcol(objs)...)
    opts = upnt(opts, "env_type", "from_file")
  elseif opts.env_type == "circle"
    isnothing(opts.separation) && error("must provide a separation between circles")
    obs = obstacles_vertical(opts.radius, opts.separation, opts.width, opts.height)
  elseif opts.env_type == "centered_circle"
    obs = obstacles_centered(opts.width, opts.height)
  elseif opts.env_type == "hex_grid"
    obs = obstacles_hexgrid(opts.radius, opts.separation, opts.width, opts.height, opts.gap)
  else
    _density = min(opts.density, 1.0)
    obs = obstacles_uniform(opts.radius, _density, opts.width, opts.height, opts.gap)
  end

  env, num = applyObstacles!(obs, opts.radius, opts.width, opts.height, min(opts.density, 1.0))
  density = sum(env .== 2) ./ length(env)
  opts = upnt(opts, "density", density)
  return env, obs[1:num], opts
end

function obstacles_hexgrid(Radius, hotspot_sep, width, height, gap=0)
  hotspots_cxy = []

  # .separation between obstacle centers is more useful
  lat_dist = hotspot_sep + 2.0 * Radius

  # .shortest dist is between horizontal neighbours
  xObstaclesNum = floor(Int, float(width) / lat_dist)
  yObstaclesNum = floor(Int, float(height) / lat_dist)

  for obsrow in range(0, yObstaclesNum)
    ypos = ceil(Int, (2 / sqrt(3)) * Radius + obsrow * lat_dist)
    for obscol in range(0, xObstaclesNum)
      xpos = ceil(Int64, Radius + (obscol + 0.5 * (obsrow % 2)) * lat_dist)
      if 1 <= xpos <= width && 1 <= ypos + gap <= height
        push!(hotspots_cxy, (xpos, ypos + gap))
      end
    end
  end
  return hotspots_cxy
end

function obstacles_uniform(radius, phi, LX, Ly, gap=0, safety=true)
  # r = 2 * radius / sqrt(3)
  r = radius

  # number density
  numberdensity = log(1 - phi) / (-π * r^2)

  ObsNum = numberdensity * (Ly - r) * LX
  ObsNum = safety ? 2 * ObsNum : ObsNum

  hotspots_cxy = []

  # .Create the list of obstacle center positions, using uniform random values
  for _ in range(1, round(Int, ObsNum))
    x = rand(1:LX)
    y = rand((gap + round(Int, r)):Ly)
    push!(hotspots_cxy, (x, y))
  end
  return hotspots_cxy
end

function obstacles_vertical(radius::Int, sep::Int, lx::Int, ly::Int)
  obstacles = NTuple{2,Int64}[]
  center = floor(Int, lx / 2)
  r = ceil(2 * radius / sqrt(3))
  foreach(xy -> push!(obstacles, xy), [(center, y) for y in (2 * r):sep:(ly - r)])
  return obstacles
end

function obstacles_centered(lx::Int, ly::Int)
  obstacles = NTuple{2,Int64}[]
  center_x = floor(Int, lx / 2)
  center_y = floor(Int, ly / 2)
  push!(obstacles, (center_x, center_y))
  return obstacles
end

function applyObstacles!(obsCenters, radius, lx, ly, area_fraction=Inf)
  env = ones(UInt8, (lx * ly))
  r = sqrt(3) * radius
  R = ceil(Int, r) + 4

  offset(x, y) = (sqrt(3) * (x - 0.5 * (y % 2)), 1 + 1.5 * (y - 1))

  _sites_covered = 0
  _obs_number = 0
  for (ox, oy) in obsCenters
    for dx in (-R):R, dy in (-R):R
      checkleq = sum(abs2, offset(ox + dx, oy + dy) .- offset(ox, oy)) <= r^2
      checkeq = isapprox(sum(abs2, offset(ox + dx, oy + dy) .- offset(ox, oy)), r^2; rtol=0.03)
      checkleq || checkeq || continue
      ny = dy + oy
      1 <= ny <= ly || continue
      nx = mod(dx + ox, 1:lx)
      env[lx * (ny - 1) + nx] == 1 && (_sites_covered += 1)
      env[lx * (ny - 1) + nx] = 2
    end

    _obs_number += 1
    _area_fraction = _sites_covered / length(env)
    if _area_fraction >= area_fraction || isapprox(_area_fraction, area_fraction; rtol=0.005)
      break
    end
  end
  return env, _obs_number
end
