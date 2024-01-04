using LinearAlgebra, Plots, DelimitedFiles, Measures
include("potentials.jl")

V = sin_two_wells

# Optimal diffusion
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

# Constant diffusion
function lp_constraint(d, I, mu_arr, p)
    return (dot(d .^ p, mu_arr .^ p) / I)^(1 / p)
end
d_constant = ones(I) / lp_constraint(ones(I), I, mu_arr, p)

# Mesh
Nmesh = 1000
XX = LinRange(0, 1, Nmesh)

# Number of tries for a given position
Ntries = 100000

# Time step
Δt = 0.01

for (d, d_string) in zip([d_opt, d_homog, d_constant], ["d_opt", "d_homog", "d_constant"])

    if !isfile(
        "data/" *
        string(V) *
        "/rejection_probability_spatial_analysis_Δt_$(Δt)_Ntries_$(Ntries)_$(d_string).txt",
    )

        println("Running for $(d_string)")
        flush(stdout)

        sqrt_2_Δt = sqrt(2 * Δt)

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

        MH_rejection_probability = Vector{Float64}(undef, Nmesh)
        for (idx, x) in enumerate(XX)
            rejection_probability = 0
            for i = 1:Ntries
                G = randn()
                d_x = d_interp(d, I, x)
                proposal = x + sqrt_2_Δt * sqrt(d_x) * G
                d_proposal = d_interp(d, I, proposal)
                V_x = V(x)
                V_proposal = V(proposal)
                sqrt_d_ratio = sqrt(d_x / d_proposal)
                G_proposal = sqrt_d_ratio * G
                alpha = (log(sqrt_d_ratio) - (V_proposal - V_x) - (G_proposal^2 - G^2) / 2)
                if log(rand()) > alpha
                    rejection_probability += 1
                end
            end
            MH_rejection_probability[idx] = rejection_probability / Ntries
        end

        writedlm(
            "data/" *
            string(V) *
            "/rejection_probability_spatial_analysis_Δt_$(Δt)_Ntries_$(Ntries)_$(d_string).txt",
            MH_rejection_probability,
        )

    else
        MH_rejection_probability = readdlm(
            "data/" *
            string(V) *
            "/rejection_probability_spatial_analysis_Δt_$(Δt)_Ntries_$(Ntries)_$(d_string).txt",
        )
    end

    if d_string == "d_opt"
        plot(XX, d_opt, label = "optimal diffusion", color = :red)
    elseif d_string == "d_homog"
        plot(XX, d_homog, label = "homogenized diffusion", color = :blue)
    elseif d_string == "d_constant"
        plot(XX, d_constant, label = "constant diffusion", color = :red)
    end
    plot!(XX, pi_arr, label = "target distribution", color = :grey)
    plot!(
        XX,
        NaN .* (XX),
        label = "MH rejection probability",
        grid = false,
        color = :green,
        right_margin = 1cm,
    )

    plot!(
        twinx(),
        XX,
        MH_rejection_probability,
        legend = false,
        xticks = :none,
        color = :green,
        y_foreground_color_axis = :green,
        y_foreground_color_border = :green,
        y_foreground_color_text = :green,
    )
    savefig(
        "data/" *
        string(V) *
        "/rejection_probability_spatial_analysis_Δt_$(Δt)_Ntries_$(Ntries)_$(d_string).png",
    )
end

# For various lower bounds
a_range = [0.2,0.4,0.6,0.8,1] #[0.05 * i for i = 0:20]
colors = palette(:balance, length(a_range))

fig = plot()
for (idx, a) in enumerate(a_range)
    local a_rounded = round(a, sigdigits = 2)
    if !isfile(
        "data/" *
        string(V) *
        "/rejection_probability_spatial_analysis_Δt_$(Δt)_Ntries_$(Ntries)_d_opt_a_$(a_rounded).txt",
    )

        println("Running for a = $(a_rounded)")
        flush(stdout)
        local dir_string = "data/" * string(V) * "/I_$(I)_a_$(a_rounded)_b_$(b_rounded)/"
        local d_opt = readdlm(dir_string * "d_opt.txt")

        sqrt_2_Δt = sqrt(2 * Δt)

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

        MH_rejection_probability = Vector{Float64}(undef, Nmesh)
        for (idx, x) in enumerate(XX)
            rejection_probability = 0
            for i = 1:Ntries
                G = randn()
                d_x = d_interp(d_opt, I, x)
                proposal = x + sqrt_2_Δt * sqrt(d_x) * G
                d_proposal = d_interp(d_opt, I, proposal)
                V_x = V(x)
                V_proposal = V(proposal)
                sqrt_d_ratio = sqrt(d_x / d_proposal)
                G_proposal = sqrt_d_ratio * G
                alpha = (log(sqrt_d_ratio) - (V_proposal - V_x) - (G_proposal^2 - G^2) / 2)
                if log(rand()) > alpha
                    rejection_probability += 1
                end
            end
            MH_rejection_probability[idx] = rejection_probability / Ntries
        end

        writedlm(
            "data/" *
            string(V) *
            "/rejection_probability_spatial_analysis_Δt_$(Δt)_Ntries_$(Ntries)_d_opt_a_$(a_rounded).txt",
            MH_rejection_probability,
        )

    else
        MH_rejection_probability = readdlm(
            "data/" *
            string(V) *
            "/rejection_probability_spatial_analysis_Δt_$(Δt)_Ntries_$(Ntries)_d_opt_a_$(a_rounded).txt",
        )
    end

    plot!(
        fig,
        XX,
        MH_rejection_probability,
        label = "a = $(a_rounded)",
        color = colors[idx],
    )
end

ylabel!("MH rejection probability")
plot!(legend = :bottomleft)
savefig(
    "data/" *
    string(V) *
    "/rejection_probability_spatial_analysis_Δt_$(Δt)_Ntries_$(Ntries)_d_opt_various_lower_bounds.png",
)
