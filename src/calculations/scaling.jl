module Observables
include("binning.jl")
include("fitting.jl")
include("../base/containers.jl")
import StatsBase: mean, var, countmap, mean_and_var, mean_and_std

function survivalFreq(graph, front, key; species=3)
  n_surviving = countmap(getfield.(graph[front], key))
  #n, v = Int64[], Int64[]

  d = Dict{String,Int64}()
  for label in 1:species
    d["n_$(label)"] = get(n_surviving, label, 0)
    d["v_$(label)"] = sum(getfield.(graph, key) .== label)
    #push!(n, get(n_surviving, label, 0))
    #push!(v, sum(getfield.(graph, key) .== label))
  end
  return namedtuple(d)
  #return (n..., v...)
end

# doc: characteristic width and heights of domains emerging from mutations
function bubbleScaling(graph, dims; delim=1, nbins=1000)
  labels = reshape(getfield.(graph, :ID_4), dims)
  maxLabels = maximum(labels)

  maxLabels > 1 || return (0.0, 0.0), Int64[]

  labelCounts = zeros(Int64, maxLabels + 1)

  bubbles_x = Array{Vector{Inf64}}(undef, length(labelCounts))
  bubbles_y = Array{Vector{Inf64}}(undef, length(labelCounts))

  # .for each height, count number sector boundaries (mutant IDs)
  for y in axes(labels, 2)
    @inbounds for x in axes(labels, 1)
      @inbounds labelCounts[1 + labels[x, y]] += 1
      push!(bubbles_x[1 + labelCounts[x, y]], y)
      push!(bubbles_y[1 + labelCounts[x, y]], y)
    end
  end

  for i in eachindex(bubbles_x, bubbles_y)
    # .skip 'small' domains
    length(bubbles_y[i]) > 10 || continue

    xmin, xmax = extrema(bubbles_x[i])
    if abs(xmin - xmax) >= dims[1] / 2
      bubbles_x[i] = map(t -> t > dims[1] / 2 ? t - L : t, bubbles_x[i])
    end

    ellipse = EllipseFitting.fitEllipse(bubbles_x[i], bubbles_y[i])
    # ismissing(ellipse) ||
  end
  #?: bin bubble sizes?
end

# doc: characteristic width and heights of domains emerging from mutations
function domainSizes(graph, dims; delim=1, nbins=1000)
  labels = reshape(getfield.(graph, :ID_4), dims)
  maxLabels = maximum(labels)

  maxLabels > 1 || return (0.0, 0.0), Int64[]

  ymax = lastindex(labels, 2)
  sectorCounts = zeros(Int64, maxLabels + 1)
  labelCounts = zeros(Int64, maxLabels + 1)

  # .for each height, count number species (mutant IDs)
  for y in axes(labels, 2)
    @inbounds for x in axes(labels, 1)
      labelCounts[1 + labels[x, y]] += 1
      if y == ymax && labels[x, y] != 0
        sectorCounts[1 + labels[x, y]] += 1
      end
    end
  end

  # largestBubble = maximum(labelCounts[sectorCounts .<= 2][2:end])

  # .find all non-zero domains, skipping the ID_4 == 0 (wild type)
  # ξ1 = filter(t -> t > 2, labelCounts[2:end]) ./ largestBubble
  ξ1 = filter(!iszero, labelCounts[2:end]) ./ reduce(*, dims)

  binmap = binMapping([0, 1]; n=nbins)
  _histogram = zeros(Int64, nbins)
  updateHist!(_histogram, ξ1, binmap)
  return mean_and_std(ξ1), _histogram
end

function sectorHeights(graph, dims)
  grid_2 = reshape(getfield.(graph, :ID_2), dims)
  sectorWidths = zeros(Int, size(grid_2))

  # .for each height, count of number of each sector (genetic ID)
  for iy in axes(grid_2, 2)
    for ix in axes(grid_2, 1)
      @inbounds sectorWidths[grid_2[ix, iy], iy] += 1
    end
  end

  ξ1 = Int[]
  ξ2 = Int[]

  for ix in axes(sectorWidths, 1)
    graph[ix].ID_3 == 1 || continue

    sector = view(sectorWidths, ix, :)

    # .ignore really small sectors bc they are susceptible to lattice effects
    peak = findfirst(iszero, sector)
    (isnothing(peak) || peak < 3) && continue

    width = maximum(sector)

    push!(ξ1, peak)
    push!(ξ2, width)
  end
  return (mean(ξ1), maximum(ξ1), mean(ξ2), maximum(ξ2))
end

function traceBoundary(graph, lx, ly, label=:ID_3; window=20, vertical=true, delim=1, stop_index=ly)
  labels = reshape(getfield.(graph, label), (lx, ly))

  # starting position, assuming split at center
  x0 = fld(lx, 2)

  vertical || (labels = transpose(labels))
  lim_h = size(labels, 1)
  lim_v = size(labels, 2)

  bndry = zeros(Float32, lim_v)

  for y in axes(labels, 2)
    y > stop_index && break
    for x in (x0 - window):(x0 + window)
      (x < 1 || x >= lim_h) && continue

      if (labels[x, y] != labels[x + 1, y]) && (labels[x, y] == delim || labels[x + 1, y] == delim)
        bndry[y] = (isodd(y) ? x - 0.5 : x)
        x0 = x
        break
      end
    end
  end
  return bndry
end

function boundaryScaling(bndry)
  variance = _roughening(bndry)
  rgn = Tuple{Int,Float64}[]
  for t in eachindex(variance)
    t0 = fld(t, 2)
    if t0 > 0 && variance[t] != 0 && variance[t0] != 0.0
      push!(rgn, (t, log(variance[t] / variance[t0]) / log(t / t0)))
    end
  end
  return rgn
end

#doc 'onloine' updates for the mean and variance of the boundary
#https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
function onlineUpdates!(existing::Matrix{<:Real}, bndry::Vector{<:Real})
  for i in eachindex(bndry)
    isfinite(bndry[i]) || continue

    (count, mean, M2) = existing[i, :]
    new_value = bndry[i]

    count += 1
    delta = new_value - mean
    mean += delta / count
    delta2 = new_value - mean
    M2 += delta * delta2
    existing[i, :] .= (count, mean, M2)
  end
  return nothing
end

#. return (sample) mean and variance from the online updates
function retrieveOnlineUpdates(existing_aggregate::Matrix{<:Real})
  count = existing_aggregate[:, 1]
  mean = existing_aggregate[:, 2]
  M2 = existing_aggregate[:, 3]
  return mean, M2 ./ count, M2 ./ (count .- 1)
end

function _roughening(bndry::Vector{<:Real})::Vector{Float64}
  variance = zeros(Float64, size(bndry))
  for t in eachindex(bndry)
    vals = bndry[1:t]
    variance[t] = var(vals)
  end
  return variance
end

function _lateralMSD(bndry::Vector{<:Real})::Vector{Float64}
  # _msd = zeros(Float64, size(bndry))
  # msd = diff(bndry .- first(bndry)) .^2
  #? determine the MSD from average of the boundary
  msd = (bndry .- first(bndry)) .^ 2
  return msd
end

function mutantFrequency(graph, width, height, label=:ID_3; target=3)
  labels = reshape(getfield.(graph, label), (width, height))
  return sum(labels .== target; dims=1)[1, :] ./ width
end

function lateralSectorSize(graph, dims)
  sectorSizes = zeros(Float64, last(dims), 3)
  grid_2 = reshape(getfield.(graph, :ID_2), dims)

  for y in axes(grid_2, 2)
    uniqueElems = countmap(view(grid_2, :, y))
    nroots = length(keys(uniqueElems))
    m, v = mean_and_var(values(uniqueElems))
    sectorSizes[y, 1] = m
    sectorSizes[y, 2] = v
    sectorSizes[y, 3] = nroots
  end
  return sectorSizes
end

end # module
