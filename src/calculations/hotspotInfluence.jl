module HotspotInfluences

function contributing!(resultContainer, graph, opts)
  r = sqrt(3) * opts.radius
  R = ceil(Int, r) + 4

  offset(x, y) = (sqrt(3) * (x - 0.5 * (y % 2)), 1 + 1.5 * (y - 1))

  for (i,(ox, oy)) in enumerate(resultContainer.htspts)
    frequency = 0
    counter = 0
    for dx in (-R):R, dy in (-R):R
      checkleq = sum(abs2, offset(ox + dx, oy + dy) .- offset(ox, oy)) <= r^2
      checkeq = isapprox(sum(abs2, offset(ox + dx, oy + dy) .- offset(ox, oy)), r^2; rtol=0.03)
      checkleq || checkeq || continue
      ny = dy + oy
      1 <= ny <= opts.height || continue
      nx = mod(dx + ox, 1:opts.width)
      if graph[opts.width * (ny - 1) + nx].ID_3 == 2
        frequency += 1 
      end
      counter += 1
    end
    resultContainer.influentialHotspots[i] += frequency > 0.5 * counter ? 1 : 0
  end
  return nothing
end

end
