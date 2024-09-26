using FileIO, Arrow, ArgParse

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
    if isfile(joinpath(root, dir, "heatmap_ID3.jld2"))
      # Get the full path of the file
      jld2_path = joinpath(root, dir, "heatmap_ID3.jld2")
      
      # Load the data from the JLD2 file
      heatmap_ID3 = load(jld2_path, "heatmap_ID3")
      
      # Save the data to Arrow format
      save_path = dirname(jld2_path)
      open(Arrow.Writer, joinpath(save_path, "heatmap_ID3.arrow")) do file
        Arrow.write(file, (heatmap_ID3=heatmap_ID3,))
      end
    end
  end
end