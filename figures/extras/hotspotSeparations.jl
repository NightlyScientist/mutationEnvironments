using CairoMakie, Distances, LinearAlgebra, StatsBase
include("../../src/base/environment.jl")
include("../calculations/OptimalPaths.jl")
include("../common/modeling.jl")

img_path = "/home/jgonzaleznunez/Projects/mutationWithLandscapes/workspace/media/images/"

λ(r, phi) = sqrt(-π * (r)^2 / log(1 - phi))
unzip(a) = map(x -> getfield.(a, x), fieldnames(eltype(a)))

L = 1000
r = 10

opts = (width=2 * L, height=L, radius=r, intensity=1, dims=(2 * L, L), ref_line=L)

sep_vs_phi = NTuple{3,Float64}[]

for phi in 0.05:0.01:0.9
  for _ in 1:5
    obs = obstacles_uniform(r, phi, L, L, 0)
    _, num = applyObstacles!(obs, r, L, L, min(phi, 1.0))
    obs = unique(obs[1:num])

    points = transpose(hcat(Float64.(first.(obs)), Float64.(last.(obs))))

    voronoiGraph = ContinuousOptPaths.voronoiTriangulation(points, opts)

    distances = Float64[]
    for (_, v) in voronoiGraph.distance
      for n in v
        push!(distances, n)
      end
    end

    #hx = getfield.(obs, 1)
    #hy = getfield.(obs, 2)
    #mat = hcat(hx, hy)

    #D = pairwise(Euclidean(), mat', mat')
    #D[diagind(D)] .= Inf
    #(distances, _) = findmin(D; dims=2)

    m, v = mean_and_std(distances)
    push!(sep_vs_phi, (phi, m, v))
  end
end

begin
  (x, y, z) = unzip(sep_vs_phi)
  fig, ax = scatter(x, y)
  ax.yscale = log10
  ax.xscale = log10
  lines!(ax, x, sqrt(2) * λ.(r, x); color=:red)
  lines!(ax, x, λ.(r, x) * 2 * sqrt(2) / sqrt(3); color=:green)
  errorbars!(ax, x, y, z; color=:black)

  display(fig)
  save(joinpath(img_path,"hotspotSeparation_$(opts.width)_$(opts.height)_$(opts.radius).png"), fig)
end
