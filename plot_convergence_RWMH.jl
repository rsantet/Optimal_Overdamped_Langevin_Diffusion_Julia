using Plots, JLD2, Statistics, LaTeXStrings
include("potentials.jl")
V = sin_two_wells
dir_string = "data/$(string(V))/convergence_RWMH/"


Nsamples = 10^4
Nsimulations = 10
Niterations = 10^6
Δt = 0.000001
nbins = 100


## Prepare plots

TV_plot = plot(yaxis=:log, xlabel="Physical time", ylabel="Total variation error")
L2_plot = plot(yaxis=:log, xlabel="Physical time", ylabel="L2 error")

d_colors = [:red, :blue, :black]
d_labels = ["optimal diffusion", "homogenized diffusion", "constant diffusion"]
XX = [i*Δt for i in 1:Niterations]

for (idx_d, d_string) in enumerate(d_labels)

    dir_string_d = dir_string * d_string * "/"

    TV_means = load_object(dir_string_d * "TV_means_$(Δt).jld2")[:]
    TV_stds = load_object(dir_string_d * "TV_stds_$(Δt).jld2")[:]
    L2_means = load_object(dir_string_d * "L2_means_$(Δt).jld2")[:]
    L2_stds = load_object(dir_string_d * "L2_stds_$(Δt).jld2")[:]

    #= plot!(
        TV_plot,
        XX, TV_means, label="$(d_labels[idx_d])", color=d_colors[idx_d]
    )
    plot!(
        L2_plot,
        XX, L2_means, label="$(d_labels[idx_d])", color=d_colors[idx_d]
    ) =#
    
    plot!(
        TV_plot,
        XX, TV_means, label="$(d_labels[idx_d])", color=d_colors[idx_d], ribbon=1.96 .* TV_stds ./ sqrt(Nsimulations)
    )
    plot!(
        L2_plot,
        XX, L2_means, label="$(d_labels[idx_d])", color=d_colors[idx_d], ribbon=1.96 .* L2_stds ./ sqrt(Nsimulations)
    )
end

# for sin_two_wells
spectral_gap_opt = 11.22782267676459966
spectral_gap_homog = 10.572299700201395
spectral_gap_constant = 0.8107051298637935


# for cos_2
#= spectral_gap_opt = 36.87891434091044118
spectral_gap_homog = 32.43337595417769
spectral_gap_constant = 30.474986632730083 =#


plot!(L2_plot,
XX, exp.(-spectral_gap_opt.*XX) / 1.75, color=:red, label=L"$\exp(-\lambda_{\mathrm{opt}}t)$", linestyle=:dash)
plot!(L2_plot,
XX, exp.(-spectral_gap_homog.*XX) / 1.75, color=:blue, label=L"$\exp(-\lambda_{\mathrm{hom}}t)$", linestyle=:dash)
plot!(L2_plot,
XX, exp.(-spectral_gap_constant.*XX) / 2, color=:black, label=L"$\exp(-\lambda_{\mathrm{cst}}t)$", linestyle=:dash, dpi=300)
ylims!(10^-1,3)
xlims!(0, 0.25)
savefig(L2_plot, dir_string * "L2_means_plot_$(Δt).png")
plot!(TV_plot,
XX, exp.(-spectral_gap_opt.*XX), color=:red, label=L"$\exp(-\lambda_{\mathrm{opt}}t)$", linestyle=:dash)
plot!(TV_plot,
XX, exp.(-spectral_gap_homog.*XX), color=:blue, label=L"$\exp(-\lambda_{\mathrm{hom}}t)$", linestyle=:dash)
plot!(TV_plot,
XX, exp.(-spectral_gap_constant.*XX), color=:black, label=L"$\exp(-\lambda_{\mathrm{cst}}t)$", linestyle=:dash)
savefig(TV_plot, dir_string * "TV_means_plot_$(Δt).png")