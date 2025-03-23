using FileIO, JLD2, DelimitedFiles

function createEnvironment!(opts)
  env_funcs = Dict(
    "circle" => obstacles_vertical,
    "centered_circle" => obstacles_centered,
    "hex_grid" => obstacles_hexgrid,
    "hot_columns" => hot_columns!,
    "uniform" => obstacles_uniform
  )

  upnt(c, k, v) = merge(c, (Symbol(k) => v,))
  opts = upnt(opts, "dims", (opts.width, opts.height))

  isnothing(opts.env_type) && error("must provide an environment type")
  haskey(env_funcs, opts.env_type) || error("invalid environment type")
  obs = env_funcs[opts.env_type](opts)

  if ~isnothing(opts.landscape)
    ispath(opts.landscape) || error("could not find environment.txt in path.")
    obs = tuple.(eachcol(readdlm(opts.landscape, ',', Int, '\n'))...)
    opts = upnt(opts, "env_type", "from_file")
  end

  env, num = addObs(obs, opts)
  opts = upnt(opts, "density", sum(env .== 2) ./ length(env))
  return env, obs[1:num], opts
end

function obstacles_hexgrid(opts)
  (; radius, separation, width, height, gap) = opts
  isnothing(separation) && error("must provide a separation between circles")

  ObsCenterList = []

  # .separation between obstacle centers is more useful
  LatticeCentreSep = separation + 2.0 * radius

  # .shortest dist is between horizontal neighbours
  xObstaclesNum = floor(Int, float(width) / LatticeCentreSep)
  yObstaclesNum = floor(Int, float(height) / LatticeCentreSep)

  for obsrow in range(0, yObstaclesNum)
    ypos = ceil(Int, (2 / sqrt(3)) * radius + obsrow * LatticeCentreSep)
    for obscol in range(0, xObstaclesNum)
      xpos = ceil(Int64, radius + (obscol + 0.5 * (obsrow % 2)) * LatticeCentreSep)
      if 1 <= xpos + gap <= width && 1 <= ypos + gap <= height
        push!(ObsCenterList, (xpos + gap, ypos + gap))
      end
    end
  end
  return ObsCenterList
end

function obstacles_uniform(opts)
  (; radius, density, width, height, gap) = opts
  # r = 2 * radius / sqrt(3)
  r = radius

  # .Units: obstacles per site. hence the 2/root(3)
  # numberdensity = log(1 - phi) / (-π * r^2 * 2 / sqrt(3))
  numberdensity = log(1 - min(density, 1.0)) / (-π * r^2)

  ObsNum = 2 * numberdensity * (height - r) * width
  ObsCenterList = []

  # .Create the list of obstacle center positions, using uniform random values
  for _ in range(1, round(Int, ObsNum))
    xy = (rand(1:width), rand((gap + round(Int, r)):height))
    push!(ObsCenterList, xy)
  end
  return ObsCenterList
end

function obstacles_vertical(opts)
  isnothing(opts.separation) && error("must provide a separation between circles")
  (; radius, separation, width, height) = opts
  obs = NTuple{2,Int64}[]
  center = floor(Int, width / 2)
  r = ceil(2 * radius / sqrt(3))
  foreach(xy -> push!(obs, xy), [(center, y) for y in (2 * r):separation:(height - r)])
  return obs
end

function obstacles_centered(opts)
  (; width, height) = opts
  obs = NTuple{2,Int64}[]
  center_x = floor(Int, width / 2)
  center_y = floor(Int, height / 2)
  push!(obs, (center_x, center_y))
  return obs
end

"""generate columns of random noise to approximate quenched noise in a 1+1D population"""
function hot_columns!(_, radius, lx, ly, phi)
  env = ones(UInt8, (lx * ly))

  # .Create the list of obstacle center positions, using uniform random values
  numberdensity = -log(1 - phi) / (radius)
  ObsNum = round(Int, numberdensity * lx)
  ObsCenterList = [rand(1:lx) for _ in 1:ObsNum]

  _obs_number = 0
  for ox in unique(ObsCenterList)
    for dx in (-radius):radius, ny in 1:ly
      nx = mod(dx + ox, 1:lx)
      env[lx * (ny - 1) + nx] = 2
    end

    _obs_number += 1
    _area_fraction = sum(env[1:lx] .== 2) / lx
    if _area_fraction >= phi || isapprox(_area_fraction, phi; rtol=0.005)
      break
    end
  end
  return env, _obs_number
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

function addObs(obs, opts)
  args = (obs, opts.radius, opts.width, opts.height, min(opts.density, 1.0))
  f = opts.env_type == "hot_columns" ? hot_columns! : applyObstacles!
  return f(args...)
end