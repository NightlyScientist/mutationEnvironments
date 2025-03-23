using FileIO, Arrow, ArgParse, JLD2
include("../common/modeling.jl")

#doc read jld file and extract lineage positions
function lineageTracks(dims, data)
  lhmap = zeros(Int64, reduce(*, dims))
  for trial in FilterTools.mask(keys(data), "trial")
    lhmap[collect(keys(data[trial]["phylogeny"]))] .+= 1
  end
  return lhmap ./ maximum(lhmap)
end

sts = ArgParseSettings()
addArgs!(sts, name; kwargs...) = add_arg_table!(sts, "--$(name)", Dict(kwargs))
addArgs!(sts, "base_path"; required=false, arg_type=String)
args = parse_args(sts)

# Get the base path from the command line argument
base_path = args["base_path"]

# Iterate over all subdirectories
for (root, dirs, files) in walkdir(base_path)
  # Check if the current directory contains the desired file
  for dir in dirs
    # Get the full path of the file
    jld2_path = joinpath(root, dir, "data_phylo.jld2")

    if isfile(jld2_path)
      dims = load(joinpath(root, dir, "Opts.jld2"), "width", "height")

      # Load the data from the JLD2 file
      file = jldopen(jld2_path)
      lineage_tracks = lineageTracks(dims, file)
      close(file)
      
      # Save the data to Arrow format
      save_path = dirname(jld2_path)
      open(Arrow.Writer, joinpath(save_path, "lineage_tracks_heatmap.arrow")) do file
        Arrow.write(file, (lineage_tracks=lineage_tracks,))
      end
    end
  end
end