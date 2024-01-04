using JLD2, LinearAlgebra, DelimitedFiles, Statistics, Plots
using Base.Threads
include("potentials.jl")

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


# Interpolation function 
function d_interp(d, I, x)
    # I is the length of d
    q = mod(x, 1)
    k = floor(Int, I * q + 1)
    k1 = k + 1
    if k == I
        k1 = 1
    elseif k == I + 1
        return d[1]
    end
    return (d[k] * (q - (k - 1) / I) + d[k1] * (k / I - q)) * I
end

# Histograms parameters
nbins = 100
nbins_XX = [i/nbins for i in 0:(nbins-1)]
gibbs_distribution = map(x->exp(-V(x)), nbins_XX)
Z_gibbs = sum(gibbs_distribution) / nbins
gibbs_distribution /= Z_gibbs

function get_idx_for_hist(q)
    tmp = mod(q, 1)
    return floor(Int64, tmp*nbins+1)
end

function normalized_hist(hist)
    return nbins*hist/sum(hist)
end

# Compute total variation distance between empirical probability distribution and Gibbs measure
function TV_distance(μ_t,μ)
    return 1/(2*nbins)*sum(abs.(μ_t-μ))
end

# Compute L^2 distance between empirical probability distribution and Gibbs measure
function L2_distance(μ_t,μ)
    return sqrt(1/nbins*sum((μ_t-μ).^2 ./ μ))
end

# Performs RWMH 
Nsamples = 10^4
Nsimulations = 10
Niterations = 10^6
Δt = 0.000001
β = 1.
x0 = 0.36544
sqrt_2_dt = sqrt(2 * Δt / β)


dir_string = "data/$(string(V))/convergence_RWMH/"
mkpath(dir_string)

d_labels = ["optimal diffusion", "homogenized diffusion", "constant diffusion"]

for (idx_d, d) in enumerate([d_opt, d_homog, d_constant])

    dir_string_d = dir_string * d_labels[idx_d] * "/"
    mkpath(dir_string_d)

    TV_values_simulations = fill(0., (Nsimulations, Niterations))
    L2_values_simulations = fill(0., (Nsimulations, Niterations))

    for idx_simulation in 1:Nsimulations

        hist_t = Array{Threads.Atomic{Float64}}(undef, (Niterations, nbins))
        for i in 1:Niterations, j in 1:nbins
            hist_t[i,j]= Threads.Atomic{Float64}(0.)
        end

        @threads for idx_sample in 1:Nsamples
            if idx_sample % 1000 == 0
                println("d = $(d_labels[idx_d]), simulation = $(idx_simulation)/$(Nsimulations), sample = $(idx_sample)/$(Nsamples)")
                flush(stdout)
            end

            Gs = randn(Niterations)
            x = rand() # Uniform distribution for the initial distribution

            idx_hist = get_idx_for_hist(x)
            Threads.atomic_add!(hist_t[1, idx_hist], 1.)

            for idx_iteration in 2:Niterations

                G = Gs[idx_iteration]
                d_x = d_interp(d, I, x)
                proposal = x + sqrt_2_dt * sqrt(d_x) * G
                d_proposal = d_interp(d, I, proposal)
                V_x = V(x)
                V_proposal = V(proposal)
                sqrt_d_ratio = sqrt(d_x / d_proposal)
                G_proposal = sqrt_d_ratio * G
                alpha = log(sqrt_d_ratio) - β*(V_proposal - V_x) - (G_proposal^2 - G^2) / 2

                if log(rand()) <= alpha
                    x = proposal
                end
                idx_hist = get_idx_for_hist(x)
                Threads.atomic_add!(hist_t[idx_iteration, idx_hist], 1.)
            end
        end

        for idx_iteration in 1:Niterations
            normed_hist = normalized_hist(
                [hist_t[idx_iteration,j][] for j in 1:nbins]
            )
            TV_values_simulations[idx_simulation, idx_iteration] = TV_distance(normed_hist, gibbs_distribution)
            L2_values_simulations[idx_simulation, idx_iteration] = L2_distance(normed_hist, gibbs_distribution)
        end
    end

    TV_means = mean(TV_values_simulations, dims=1)
    L2_means = mean(L2_values_simulations, dims=1)
    TV_stds = std(TV_values_simulations, dims=1)
    L2_stds = std(L2_values_simulations, dims=1)

    save_object(dir_string_d * "TV_means_$(Δt).jld2", TV_means)
    save_object(dir_string_d * "TV_stds_$(Δt).jld2", TV_stds)
    save_object(dir_string_d * "L2_means_$(Δt).jld2", L2_means)
    save_object(dir_string_d * "L2_stds_$(Δt).jld2", L2_stds)
end

