using CairoMakie, Distributions, StatsBase 
include("../../src/calculations/fitting.jl")

alpha = 5
beta = 3
N = 500
DIM = 2

StatsBase.Random.seed!(2)

# .Generate random points on the unit circle by sampling uniform angles
theta = 2 * π * rand(Float64, (N, 1))
eps_noise = 0.2 * rand(Normal(), (N, 1))
circle = hcat(cos.(theta), sin.(theta))

# .Stretch and rotate circle to an ellipse with random linear tranformation
B = rand(-3:3, (DIM, DIM))
noisy_ellipse = circle * B .+ eps_noise

# .Extract x coords and y coords of the ellipse as column vectors
x = noisy_ellipse[:, 1] .+ 10
y = noisy_ellipse[:, 2] .+ 10

ellipse = EllipseFitting.fitEllipse(x,y)
ex, ey = EllipseFitting.get_ellipse_pts(ellipse)

# .Plot the noisy data
fig, ax = scatter(X, Y; label="Data Points")

# .Plot the original ellipse from which the data was generated
phi = LinRange(0, 2π, 1000)
c = hcat(cos.(phi), sin.(phi))
ground_truth_ellipse = c * B .+ 10

scatter!(ax, ground_truth_ellipse[:, 1], ground_truth_ellipse[:, 2]; color=:black, label="generating ellipse")

ex, ey = get_ellipse_pts(params; npts=1000)
scatter!(ax, ex, ey; color=:red)

display(fig)