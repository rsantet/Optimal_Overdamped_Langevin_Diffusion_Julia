using DelimitedFiles, Plots
include("potentials.jl")
include("main.jl")

# Parameters for the optimization algorithm
I = 1000
p = 2
b = Inf
V = sin_two_wells

function mu(x)
    return exp(-V(x))
end

XX = [i / I for i = 0:I-1]
mu_arr = map(x -> exp(-V(x)), XX)
Z = sum(mu_arr) / I
pi_arr = mu_arr / Z
d_homog = 1 ./ mu_arr

# Preparing plots
fig = plot()
fig_normalized = plot()
gaps = []

a_range = [0, 0.2, 0.4, 0.6, 0.8, 1] #[0.05 * i for i = 0:20]
colors = palette(:balance, length(a_range))
for (idx, a) in enumerate(a_range)

    a_rounded = round(a, sigdigits = 2)
    b_rounded = round(b, sigdigits = 2)
    dir_string = "data/" * string(V) * "/I_$(I)_a_$(a_rounded)_b_$(b_rounded)/"

    if !isfile(dir_string * "first_eigenvalue.txt")
        optim_algo(V, I, mu_arr, p = p, a = a, b = b)
    end

    d_opt = readdlm(dir_string * "d_opt.txt")
    plot!(fig, XX, d_opt, label = "a=$(a_rounded)", color = colors[idx])
    plot!(
        fig_normalized,
        XX,
        d_opt ./ d_homog,
        label = "a=$(a_rounded)",
        color = colors[idx],
    )

    d_opt_gap = readdlm(dir_string * "d_opt_gap.txt")[1]
    println(d_opt_gap)
    push!(gaps, d_opt_gap)
end
plot!(fig, legend = :topright)
plot!(fig_normalized, legend = :bottomleft)

savefig(fig, "data/" * string(V) * "/varying_lower_bound_optimal_diffusions.png")
savefig(
    fig_normalized,
    "data/" * string(V) * "/varying_lower_bound_optimal_diffusions_normalized.png",
)


plot(a_range, gaps, label = "", yaxis = :log)
xlabel!("a")
ylabel!("Spectral gaps")
savefig("data/" * string(V) * "/varying_lower_bound_spectral_gaps.png")
