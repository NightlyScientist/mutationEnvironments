module SimulationSettings
using Parameters: Parameters
using ArgParse: ArgParse
include("../base/containers.jl")

Parameters.@kwdef struct Setting
  name::Union{String,Missing} = missing
  required::Bool = false
  #arg_type::Union{DataType,Missing} = missing
  arg_type::Union{DataType,Missing} = missing
  default::Any = missing
  abbr::Union{String,Missing} = missing
end

Parameters.@kwdef mutable struct Settings
  program_name::Union{String,Missing} = missing
  options::Union{Vector{Setting},Missing} = missing
  exclusive::Union{Vector{NTuple{2,String}},Missing} = missing
  flags::Union{Vector{NTuple{2,String}},Missing} = missing
end

ArgParse.parse_item(::Type{NTuple{3,T}}, x::AbstractString) where {T} =
  Tuple(convert.(T, parse.(Float64, split(x, ','))))

"""parse settings from command line using the specification in simmulation_settings"""
function parse_settings(settings::Settings)
  _add!(sts, name; kwargs...) = ArgParse.add_arg_table!(sts, "--$(name)", Dict(kwargs))

  sts = ArgParse.ArgParseSettings()

  # model variable that isn't used after including model.jl
  _add!(sts, "model"; required=false, arg_type=String)

  # general settings
  for s in settings.options
    # check if default value provide if not required
    # refactor: require default value for all non-required values
    # (ismissing(s.default) && s.required) && (error("non-required settings needs a default value"))

    if ismissing(s.default)
      _add!(sts, s.name; required=s.required, arg_type=s.arg_type)
    else
      _add!(sts, s.name; default=s.default, required=s.required, arg_type=s.arg_type)
    end
  end

  # flags with :store_true
  foreach(flag -> _add!(sts, flag; action=:store_true), getfield.(settings.flags, 1))

  # exclusive option groups
  ArgParse.add_arg_group!(sts, "exlusive"; exclusive=true)
  foreach(excl -> _add!(sts, excl; action=:store_true), getfield.(settings.exclusive, 1))
  return namedtuple(ArgParse.parse_args(sts))
end

"""get all settings that have an abbreviation"""
function get_abbrs(settings::Settings, parsed::Union{Missing,NamedTuple}=missing)
  _d = Dict(getfield.(settings.options, :name) .=> getfield.(settings.options, :abbr))
  _d = merge(_d, get_abbrs(settings.exclusive), get_abbrs(settings.flags))
  return filter_optionals(_d, parsed)
end

get_abbrs(settings::Vector{NTuple{2,String}}) = Dict(getfield.(settings, 1) .=> getfield.(settings, 2))

function filter_optionals(options::Dict, parsed::Union{Missing,NamedTuple})
  ismissing(parsed) && return filter(kv -> ~ismissing(last(kv)) && last(last(kv)) != "", options)
  fltr(k) = (haskey(parsed, k) && typeof(getfield(parsed, k)) <: Bool && return getfield(parsed, k)) || return true
  return filter(kv -> ~ismissing(last(kv)) && last(kv) != "" && fltr(Symbol(first(kv))), options)
end

function print_settings(parsedArgs::NamedTuple)
  println("Simulation Info:")
  foreach(k -> println("  > $k  =>  $(parsedArgs[k])"), keys(parsedArgs))
end

function ensurePath(cli::Settings, parsed::NamedTuple)::NamedTuple
  abbrs = get_abbrs(cli, parsed)
  opts = collect(keys(abbrs))[sortperm(collect(values(abbrs)); by=first)]

  parsedOpts = []
  appendOpts = []
  for opt in opts
    haskey(parsed, Symbol(opt)) || continue
    val = getfield(parsed, Symbol(opt))
    eltype(val) <: Float64 && (val = round.(val, digits=3))
    if eltype(val) <: Bool
      push!(appendOpts, "$(abbrs[opt])")
    elseif ~isnothing(val)
      push!(parsedOpts, "$(abbrs[opt])_$(val)")
    end
  end

  parsedOpts = append!(parsedOpts, appendOpts)

  program_name = ismissing(cli.program_name) ? "" : "$(cli.program_name)_"
  path = mkpath(parsed.outputPath * "/$(program_name)" * join(parsedOpts, ","))

  if ispath(path) && ~isempty(readdir(path))
    ~parsed.rewrite && (@info " ! output path not empty"; exit())
    for ext in [".txt", ".jld2", ".arrow", ".csv", ".png", ".md"]
      foreach(rm, filter(endswith(ext), readdir(path; join=true)))
    end
  end

  # write ARGS to log file
  log_dir = joinpath(path, "logs")
  isdir(log_dir) || mkdir(log_dir)
  log_file = joinpath(log_dir, "log.md")

  open(log_file, "w") do file
    write(file, "command line: $(join(ARGS, ' '))\n\n")
    write(file, "Simulation Info:\n")
    foreach(k -> write(file, "  > $k  =>  $(parsed[k])\n"), keys(parsed))
  end
  return merge(parsed, (outputPath=path,))
end

end
