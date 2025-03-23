function setup_active(; num=4)
  active = Vector{HashMaps.HashVec{UInt32}}(undef, num)
  foreach(i -> active[i] = HashMaps.HashVec{UInt32}(), collect(1:1:num))
  return active
end

function linear_initialization(opts)
  haskey(opts, :standing_variation) || error("opts does not contain key :standing_variation")
  haskey(opts, :n_species) || error("opts does not contain key :n_species")

  strainID = ones(Int, opts.dims[1])

  if opts.standing_variation
    # default is to have a uniform distribution
    strainID = rand([1, opts.n_species], opts.dims[1])

    if contains(opts.initial_type, "split")
      strainID = ones(Int, opts.dims[1])
      strainID[(cld(opts.dims[1], 2) + 1):end] .= opts.n_species
    elseif contains(opts.initial_type, "alt")
      types = [1, opts.n_species]
      strainID = [types[mod(x, 1:2)] for x in 1:opts.dims[1]]
    end
  end
  return strainID
end

function initializePopulation!(graph, cntns, env, opts, container, save_state)
  haskey(opts, :fixed_initializations) || error("opts does not contain key :fixed_initializations")
  populate_linear!(graph, cntns, env, opts, container, save_state)
end

function populate_linear!(graph, cntns, env, opts, container, save_state; row=1)
  haskey(opts, :dims) || error("opts does not contain key :dims")

  active = setup_active(; num=opts.n_species * 2)

  initial_ids = linear_initialization(opts)

  # keep previous config if option is set
  if opts.fixed_initializations
    if ismissing(save_state.initializations)
      save_state.initializations = initial_ids
    else
      initial_ids = save_state.initializations
    end
  end
    
  for col in 1:opts.dims[1]
    #. set group affliliation with strain type (f, s, b), matching the environment
    nodeID = initial_ids[col]

    # shift group affliation to match environment
    env[col] == 2 && (nodeID += opts.n_species)

    nodeIndx = opts.dims[1] * (row - 1) + col
    graph[nodeIndx].filled = true
    graph[nodeIndx].ID_1 = nodeID
    graph[nodeIndx].ID_2 = col
    graph[nodeIndx].ID_3 = initial_ids[col]

    # don't add site to front if it already has no empty nbors
    if graph[nodeIndx].nbors > 0
      HashMaps.add!(active[nodeID], nodeIndx)
      container[nodeID] += 1
    end

    updateNeighbors!(graph, cntns, nodeIndx, active, container)
  end
  return active
end