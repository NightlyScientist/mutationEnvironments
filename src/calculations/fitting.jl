module EllipseFitting
using LinearAlgebra

#> https://scipython.com/blog/direct-linear-least-squares-fitting-of-an-ellipse/
# .implementation of the direct least squares method for fitting an ellipse to a set of points in the plane

struct Ellipse
  x0::Float64
  y0::Float64
  ap::Float64
  bp::Float64
  e::Float64
  phi::Float64
end

function solution(x, y)::Vector{Float64}
  D1 = hcat(x .^ 2, x .* y, y .^ 2)
  D2 = hcat(x, y, ones(length(x)))
  S1 = D1' * D1
  S2 = D1' * D2
  S3 = D2' * D2
  T = -1 .* inv(S3) * S2'
  M = S1 + S2 * T
  C = [0 0 2; 0 -1 0; 2 0 0]
  M = inv(C) * M

  # eigval = eigvals(M)
  eigvec = eigvecs(M)
  con = @. 4 * eigvec[:, 1] * eigvec[:, 3] - eigvec[:, 2]^2
  ak = eigvec[:, findall(!iszero, con .> 0)[1]]
  return vcat(ak, T * ak)
end

function canonical(coeffs)::Union{Missing,Ellipse}
  a = coeffs[1]
  b = coeffs[2] / 2
  c = coeffs[3]
  d = coeffs[4] / 2
  f = coeffs[5] / 2
  g = coeffs[6]

  den = b^2 - a * c
  if den > 0
    @warn "coeffs do not represent an ellipse: b^2 - 4ac < 0"
    return missing
  end

  # The location of the ellipse centre.
  x0, y0 = (c * d - b * f) / den, (a * f - b * d) / den

  num = 2 * (a * f^2 + c * d^2 + g * b^2 - 2 * b * d * f - a * c * g)
  fac = sqrt((a - c)^2 + 4 * b^2)

  # The semi-major and semi-minor axis lengths (these are not sorted).
  ap = sqrt(num / den / (fac - a - c))
  bp = sqrt(num / den / (-fac - a - c))

  # Sort the semi-major and semi-minor axis lengths but keep track of
  # the original relative magnitudes of width and height.
  width_gt_height = true
  if ap < bp
    width_gt_height = false
    ap, bp = bp, ap
  end

  # The eccentricity.
  r = (bp / ap)^2
  if r > 1
    r = 1 / r
  end
  e = sqrt(1 - r)

  # The angle of anticlockwise rotation of the major-axis from x-axis.
  if b == 0
    phi = a < c ? 0 : π / 2
  else
    phi = atan((2 * b) / (a - c)) / 2
    if a > c
      phi += π / 2
    end
  end

  if !width_gt_height
    # Ensure that phi is the angle to rotate to the semi-major axis.
    phi += π / 2
  end
  phi = phi % π

  return Ellipse(x0, y0, ap, bp, e, phi)
end

function get_ellipse_pts(ellipse::Ellipse; npts=100, tmin=0, tmax=2 * π)
  t = LinRange(tmin, tmax, npts) 

  x = @. ellipse.x0 + ellipse.ap * cos(t) * cos(ellipse.phi) - ellipse.bp * sin(t) * sin(ellipse.phi)
  y = @. ellipse.y0 + ellipse.ap * cos(t) * sin(ellipse.phi) + ellipse.bp * sin(t) * cos(ellipse.phi)
  return x, y
end

function fitEllipse(x, y)::Union{Missing,Ellipse}
  sol = solution(x, y)
  return canonical(sol)
end

end