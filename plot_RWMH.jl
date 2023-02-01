using DelimitedFiles, Plots, LinearAlgebra
include("potentials.jl")
include("RWMH.jl")

# Specifying the potential
V = sin_two_wells

# Get optimal diffusion from optimization algorithm
I = 500
p = 2
a = 0.0
b = Inf
a_rounded = round(a, sigdigits=2)
b_rounded = round(b, sigdigits=2)
dir_string = "data/" * string(V) * "/I_$(I)_a_$(a_rounded)_b_$(b_rounded)/"
d_opt = readdlm(dir_string * "d_opt.txt")

# Set optimal diffusion in the homogenized limit
XX = [i / I for i in 0:I-1]
mu_arr = map(x -> exp(-V(x)), XX)
Z = sum(mu_arr) / I
pi_arr = mu_arr / Z

d_homog = 1.0 ./ mu_arr

"""
    lp_constraint(d, I, mu_arr, p)

Compute the L^p constraint for the diffusion coefficient d

# Arguments
- `d::Vector`: the diffusion coefficient values for the point in the mesh
- `I::Int`: number of point in the mesh (length of vector d)
- `mu_arr::Vector`: the approximation of the unnormalized density of the Gibbs measure
- `p::Real`: L^p constraint.
"""
function lp_constraint(d, I, mu_arr, p)
    return (dot(d .^ p, mu_arr .^ p) / I)^(1 / p)
end
d_constant = ones(I) / lp_constraint(ones(I), I, mu_arr, p)

# Performs RWMH simulations using the same Gaussian increments each time
Δt = 0.0001
x0 = 0.
N_it = 10^6
Gs = randn(N_it)
β = 1.

println("Parameters:\n Δt=$(Δt),\nx0=$(x0),\nN_it=$(N_it),\nβ=$(β)")
println("RWMH for optimal diffusion coefficient")
trajectory_opt, MH_opt = RWMH(
    d_opt,
    I, Δt, N_it, x0, V, β, Gs
)

println("RWMH for homogenized diffusion coefficient")
trajectory_homog, MH_homog = RWMH(
    d_homog,
    I, Δt, N_it, x0, V, β, Gs
)

println("RWMH for constant diffusion coefficient")
trajectory_constant, MH_constant = RWMH(
    d_constant,
    I, Δt, N_it, x0, V, β, Gs
)

println("MH rejection probabilies:")
println(
    "Optimal: $(MH_opt)\nHomogenized: $(MH_homog)\nConstant: $(MH_constant)"
)
# MH rejection probabilies for Δt = 0.01, V = sin_two_wells
# Optimal: 0.312776
# Homogenized: 0.279466
# Constant: 0.324239

# MH rejection probabilies for Δt = 0.001, V = sin_two_wells
# Optimal: 0.172784
# Homogenized: 0.117015
# Constant: 0.114145

# MH rejection probabilies for Δt = 0.0001, V = sin_two_wells
# Optimal: 0.064226
# Homogenized: 0.040612
# Constant: 0.037898

# Plot each trajectory
plot(
    1:N_it,
    trajectory_opt,
    label="optimal diffusion",
    color=:red
)
plot!(
    1:N_it,
    trajectory_homog,
    label="homogenized diffusion",
    color=:blue
)
plot!(
    1:N_it,
    trajectory_constant,
    label="constant diffusion",
    color=:black
)
savefig(dir_string * "trajectories.png")

# Histograms
h_opt = histogram(
    mod.(trajectory_opt,1),
    label="optimal diffusion",
    color=:red,
    normalize=:pdf
)
plot!(
    XX, pi_arr,
    label="Target distribution",
    linewidth=4,
    color=:grey
)
h_homog = histogram(
    mod.(trajectory_homog,1),
    label="homogenized diffusion",
    color=:blue,
    normalize=:pdf
)
plot!(
    XX, pi_arr,
    label="Target distribution",
    linewidth=4,
    color=:grey
)
h_constant = histogram(
    mod.(trajectory_constant,1),
    label="constant diffusion",
    color=:black,
    normalize=:pdf
)
plot!(
    XX, pi_arr,
    label="Target distribution",
    linewidth=4,
    color=:grey
)
plot(
    h_opt, h_homog, h_constant, layout = (3,1)
)
savefig(dir_string * "histograms.png")