using Plots, DelimitedFiles, LinearAlgebra, Statistics
include("potentials.jl")
include("RWMH.jl")


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

# MSD simulation parameters
N_sim = 10^3
N_it = 10^6
Δt = 10^(-4)
x0 = 0
β = 1

# Plot preparation
linestyles = [:solid, :dash, :dashdot]
linewidths = [3, 2, 2]
colors = [:red, :blue, :black]
fig = plot()

# Specifying the potential
V = sin_two_wells

# Get optimal diffusion from optimization algorithm
I = 1000
p = 2
a = 0.0
b = Inf
a_rounded = round(a, sigdigits = 2)
b_rounded = round(b, sigdigits = 2)
dir_string = "data/" * string(V) * "/I_$(I)_a_$(a_rounded)_b_$(b_rounded)/"
d_opt = readdlm(dir_string * "d_opt.txt")

# Set optimal diffusion in the homogenized limit
XX = [i / I for i = 0:I-1]
mu_arr = map(x -> exp(-V(x)), XX)
Z = sum(mu_arr) / I
d_homog = 1.0 ./ mu_arr

# constant diffusion
d_constant = ones(I) / lp_constraint(ones(I), I, mu_arr, p)

for (idx, (d, d_string)) in
    enumerate(zip([d_opt, d_homog, d_constant], ["d_opt", "d_homog", "d_constant"]))

    if !isfile("data/" * string(V) * "/msd_means_$(d_string).txt")
        msd = zeros((N_sim, N_it))

        for n = 1:N_sim
            println("Running msd iteration $(n)/$(N_sim) for $(d_string)")
            flush(stdout)
            trajectory, _ = RWMH(d, I, Δt, N_it, x0, V, β)
            msd[n, :] = cumsum(trajectory .^ 2) ./ collect(1:N_it)
        end

        msd_means = mean(msd, dims = 1)[1, :]
        msd_stds = std(msd, dims = 1)[1, :]

        writedlm("data/" * string(V) * "/msd_means_$(d_string).txt", msd_means)
        writedlm("data/" * string(V) * "/msd_stds_$(d_string).txt", msd_stds)

    else
        msd_means = readdlm("data/" * string(V) * "/msd_means_$(d_string).txt")
        msd_stds = readdlm("data/" * string(V) * "/msd_stds_$(d_string).txt")
    end
    plot!(
        fig,
        1:N_it,
        msd_means,
        label = d_string,
        linestyle = linestyles[idx],
        linewidth = linewidths[idx],
        color = colors[idx],
    )
end

xlabel!("Number of iterations")
ylabel!("Mean squared distance")
savefig("data/" * string(V) * "/msd.png")


# For various lower bounds
a_range = [0, 0.2, 0.4, 0.6, 0.8, 1] #[0.05 * i for i = 0:20]
colors = palette(:balance, length(a_range))

fig = plot()
for (idx, a) in enumerate(a_range)
    local a_rounded = round(a, sigdigits = 2)
    if !isfile("data/" * string(V) * "/msd_means_d_opt_a_$(a_rounded).txt")
        local dir_string = "data/" * string(V) * "/I_$(I)_a_$(a_rounded)_b_$(b_rounded)/"
        local d_opt = readdlm(dir_string * "d_opt.txt")

        msd = zeros((N_sim, N_it))

        for n = 1:N_sim
            println("Running msd iteration $(n)/$(N_sim) for a = $(a_rounded)")
            flush(stdout)
            trajectory, _ = RWMH(d_opt, I, Δt, N_it, x0, V, β)
            msd[n, :] = cumsum(trajectory .^ 2) ./ collect(1:N_it)
        end

        msd_means = mean(msd, dims = 1)[1, :]
        msd_stds = std(msd, dims = 1)[1, :]

        writedlm("data/" * string(V) * "/msd_means_d_opt_a_$(a_rounded).txt", msd_means)
        writedlm("data/" * string(V) * "/msd_stds_d_opt_a_$(a_rounded).txt", msd_stds)

    else
        msd_means = readdlm("data/" * string(V) * "/msd_means_d_opt_a_$(a_rounded).txt")
        msd_stds = readdlm("data/" * string(V) * "/msd_stds_d_opt_a_$(a_rounded).txt")
    end

    plot!(fig, 1:N_it, msd_means, label = "a = $(a_rounded)", color = colors[idx])
end

xlabel!("Number of iterations")
ylabel!("Mean squared distance")
plot!(legend=:topleft)
savefig("data/" * string(V) * "/msd_various_lower_bounds.png")