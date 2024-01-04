using Plots, DelimitedFiles, LinearAlgebra, Statistics
include("potentials.jl")

Δt_range = [0.00005, 0.0001, 0.0005, 0.001, 0.005, 0.01]
linestyles = [:solid, :dash, :dashdot, :dashdotdot]
linewidths = [3, 2, 2]
colors = [:red, :blue, :black]
labels = ["optimal diffusion", "homogenized diffusion", "constant diffusion"]
V = sin_two_wells

N = 10000

T = plot()
MH = plot()

for (idx, d_string) in enumerate(["d_opt", "d_homog", "d_constant"])
    if !isfile("data/$(string(V))/means_transition_times_N_$(N)_$(d_string).txt")

        T_means = []
        MH_means = []
        T_stds = []
        MH_stds = []

        for Δt in Δt_range
            file_string =
                "data/" * string(V) * "/transition_times_Δt_$(Δt)_N_$(N)_$(d_string)"
            times = readdlm(file_string * ".txt")
            MH_ratios = readdlm(file_string * "_MH_ratios.txt")
            T_mean = mean(times)
            MH_mean = mean(MH_ratios)
            T_std = std(times, mean = T_mean)
            MH_std = std(MH_ratios, mean = MH_mean)
            push!(T_means, T_mean)
            push!(T_stds, T_std)
            push!(MH_means, MH_mean)
            push!(MH_stds, MH_std)
        end

        writedlm(
            "data/" * string(V) * "/means_transition_times_N_$(N)_$(d_string).txt",
            T_means,
        )
        writedlm(
            "data/" * string(V) * "/means_MH_transition_times_N_$(N)_$(d_string).txt",
            MH_means,
        )
        writedlm(
            "data/" * string(V) * "/stds_transition_times_N_$(N)_$(d_string).txt",
            T_stds,
        )
        writedlm(
            "data/" * string(V) * "/stds_MH_transition_times_N_$(N)_$(d_string).txt",
            MH_stds,
        )

    else
        T_means =
            readdlm("data/" * string(V) * "/means_transition_times_N_$(N)_$(d_string).txt")
        MH_means = readdlm(
            "data/" * string(V) * "/means_MH_transition_times_N_$(N)_$(d_string).txt",
        )
        T_stds =
            readdlm("data/" * string(V) * "/stds_transition_times_N_$(N)_$(d_string).txt")
        MH_stds = readdlm(
            "data/" * string(V) * "/stds_MH_transition_times_N_$(N)_$(d_string).txt",
        )
    end

    plot!(
        T,
        Δt_range,
        T_means,
        yerror = T_stds .* 1.96 ./ sqrt(N),
        label = labels[idx],
        color = colors[idx],
        linewidth = linewidths[idx],
        linestyle = linestyles[idx],
        xaxis = :log,
        legend = :topleft,
    )
    xlabel!("Δt")
    ylabel!("Physical time")

    plot!(
        MH,
        Δt_range,
        MH_means,
        yerror = MH_stds .* 1.96 ./ sqrt(N),
        label = labels[idx],
        color = colors[idx],
        linewidth = linewidths[idx],
        linestyle = linestyles[idx],
        xaxis = :log,
        yaxis = :log,
        legend = :bottomright,
    )
end

plot!(MH, Δt_range, 10 .* Δt_range .^ (1 / 2), label = "Δt^{1/2}")
xlims!(0.00005, 0.01)
ylims!(10^(-2), 1)
xlabel!("Δt")
ylabel!("Rejection probability")

savefig(T, "data/" * string(V) * "/transition_times_N_$(N).png")
savefig(MH, "data/" * string(V) * "/transition_times_MH_ratios_N_$(N).png")

# For various lower bounds
a_range = [0, 0.2, 0.4, 0.6, 0.8, 1] #[0.05 * i for i = 0:20]
colors = palette(:curl, length(a_range))
T = plot()
MH = plot()

for (idx, a) in enumerate(a_range)
    local a_rounded = round(a, sigdigits = 2)

    if !isfile("data/$(string(V))/means_transition_times_N_$(N)_d_opt_a_$(a_rounded).txt")

        T_means = []
        MH_means = []
        T_stds = []
        MH_stds = []

        for Δt in Δt_range
            file_string =
                "data/" *
                string(V) *
                "/transition_times_Δt_$(Δt)_N_$(N)_d_opt_a_$(a_rounded)"
            times = readdlm(file_string * ".txt")
            MH_ratios = readdlm(file_string * "_MH_ratios.txt")
            T_mean = mean(times)
            MH_mean = mean(MH_ratios)
            T_std = std(times, mean = T_mean)
            MH_std = std(MH_ratios, mean = MH_mean)
            push!(T_means, T_mean)
            push!(T_stds, T_std)
            push!(MH_means, MH_mean)
            push!(MH_stds, MH_std)
        end

        writedlm(
            "data/" * string(V) * "/means_transition_times_N_$(N)_d_opt_a_$(a_rounded).txt",
            T_means,
        )
        writedlm(
            "data/" *
            string(V) *
            "/means_MH_transition_times_N_$(N)_d_opt_a_$(a_rounded).txt",
            MH_means,
        )
        writedlm(
            "data/" * string(V) * "/stds_transition_times_N_$(N)_d_opt_a_$(a_rounded).txt",
            T_stds,
        )
        writedlm(
            "data/" *
            string(V) *
            "/stds_MH_transition_times_N_$(N)_d_opt_a_$(a_rounded).txt",
            MH_stds,
        )

    else
        T_means = readdlm(
            "data/" * string(V) * "/means_transition_times_N_$(N)_d_opt_a_$(a_rounded).txt",
        )
        MH_means = readdlm(
            "data/" *
            string(V) *
            "/means_MH_transition_times_N_$(N)_d_opt_a_$(a_rounded).txt",
        )
        T_stds = readdlm(
            "data/" * string(V) * "/stds_transition_times_N_$(N)_d_opt_a_$(a_rounded).txt",
        )
        MH_stds = readdlm(
            "data/" *
            string(V) *
            "/stds_MH_transition_times_N_$(N)_d_opt_a_$(a_rounded).txt",
        )
    end

    plot!(
        T,
        Δt_range,
        T_means,
        yerror = T_stds .* 1.96 ./ sqrt(N),
        label = "optimal diffusion, a=$(a_rounded)",
        xaxis = :log,
        legend = :topleft,
        color = colors[idx],
    )
    xlabel!("Δt")
    ylabel!("Physical time")

    plot!(
        MH,
        Δt_range,
        MH_means,
        yerror = MH_stds .* 1.96 ./ sqrt(N),
        label = "optimal diffusion, a=$(a_rounded)",
        xaxis = :log,
        yaxis = :log,
        legend = :bottomright,
        color = colors[idx],
    )
end

plot!(MH, Δt_range, 10 .* Δt_range .^ (1 / 2), label = "Δt^{1/2}")
xlims!(0.00005, 0.01)
ylims!(10^(-2), 1)
xlabel!("Δt")
ylabel!("Rejection probability")

savefig(T, "data/" * string(V) * "/transition_times_N_$(N)_various_lower_bounds.png")
savefig(
    MH,
    "data/" * string(V) * "/transition_times_MH_ratios_N_$(N)_various_lower_bounds.png",
)
