using CairoMakie, FileIO, JLD2
include("../common/modeling.jl")

base_path = "/home/jgonzaleznunez/Projects/mutationWithLandscapes/workspace/simulations/speed_relation"

let
  fig = Figure()
  ax = Axis(fig[1, 1])

  for (root, dirs, files) in walkdir(base_path)
    for dir in dirs
      extra_data = joinpath(root, dir, "data_extras.jld2")
      opts_file = joinpath(root, dir, "inputs.jld2")
      isfile(extra_data) || continue

      opts = load(opts_file, "cli")
      println(opts.width, " ", opts.height)
      jldopen(extra_data) do file
        for k in keys(file)
          occursin("trial", k) || continue
          t = getfield.(file[k]["front"], 1)
          m = getfield.(file[k]["front"], 2)

          fit = FitModels.linearfit(t, m)
          scatter!(ax, opts.selection, fit.a, color="black")
        end
      end
    end
  end

  x = LinRange(0, 0.1, 100)
  lines!(ax, x, 1.578 .- (sqrt(2) * 1.154)*x, linewidth=10, color="green")
  display(fig)
end