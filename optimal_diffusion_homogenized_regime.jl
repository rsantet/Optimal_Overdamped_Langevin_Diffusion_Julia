using Plots
include("potentials.jl")
include("main.jl")


function get_optimal_diffusion_k(
    V,
    I,
    k;
    p=2.0,
    a=0.0,
    b=Inf,
    tol=10^(-1),
    max_it=200
)

    function Vk(q)
        return V(k * q)
    end

    XX = [i / I for i in 0:(I-1)]
    mu_arr = map(x -> exp(-Vk(x)), XX)
    Z = sum(mu_arr) / I
    pi_arr = mu_arr / Z

    d_opt, d_opt_gap = optim_algo(
        Vk,
        I,
        mu_arr;
        p,
        a,
        b,
        tol=tol,
        max_it=max_it,
        save=false,
        rewrite_save=false
    )

    a_rounded = round(a, sigdigits=2)
    b_rounded = round(b, sigdigits=2)
    dir_string = "data/" * string(V) * "/I_$(I)_a_$(a_rounded)_b_$(b_rounded)/optimal_diffusion_homogenized_regime_k_$(k)/"
    mkpath(dir_string)
    d_opt_min = minimum(d_opt)

    mkpath(
        dir_string
    )

    writedlm(
        dir_string * "d_opt.txt", d_opt
    )
    writedlm(
        dir_string * "d_opt_gap.txt", [d_opt_gap]
    )
    writedlm(
        dir_string * "d_opt_min.txt", [d_opt_min]
    )
    return d_opt, d_opt_gap, d_opt_min
end


function optimal_diffusion_homogenized_regime(
    V, I, k_range, a, b
)

    # Performs computation if necessary
    for k in k_range
        a_rounded = round(a, sigdigits=2)
        b_rounded = round(b, sigdigits=2)
        dir_string = "data/" * string(V) * "/I_$(I)_a_$(a_rounded)_b_$(b_rounded)/optimal_diffusion_homogenized_regime_k_$(k)/"
        if !isfile(dir_string * "d_opt.txt")
            d_opt, d_opt_gap, d_opt_min = get_optimal_diffusion_k(
                V, I * k, k
            )
            mkpath(dir_string)
            writedlm(dir_string * "d_opt.txt", d_opt)
            writedlm(dir_string * "d_opt_gap.txt", d_opt_gap)
            writedlm(dir_string * "d_opt_min.txt", d_opt_min)
        end
    end

    XX = [i / I for i in 0:(I-1)]
    mu_arr = map(x -> exp(-V(x)), XX)
    Z = sum(mu_arr) / I
    d_homog = 1 ./ mu_arr

    # Plot diffusion coefficient on rescaled unit cell
    plot()
    for k in k_range
        a_rounded = round(a, sigdigits=2)
        b_rounded = round(b, sigdigits=2)
        dir_string = "data/" * string(V) * "/I_$(I)_a_$(a_rounded)_b_$(b_rounded)/optimal_diffusion_homogenized_regime_k_$(k)/"
        d_opt = readdlm(dir_string * "d_opt.txt")
        plot!(
            XX, d_opt[:][1:I], label="k = $(k)"
        )
    end
    plot!(
        XX,
        d_homog,
        label="homogenized diffusion",
        color=:blue
    )
    savefig("data/$(string(V))/optimal_diffusion_homogenized_regime_I_$(I)_d_opt.png")

    # Plot L^2 error between optimal diffusion and homogenized diffusion
    plot()
    for k in k_range
        a_rounded = round(a, sigdigits=2)
        b_rounded = round(b, sigdigits=2)
        dir_string = "data/" * string(V) * "/I_$(I)_a_$(a_rounded)_b_$(b_rounded)/optimal_diffusion_homogenized_regime_k_$(k)/"
        d_opt = readdlm(dir_string * "d_opt.txt")
        
        L2_error = sqrt(sum((d_homog - d_opt[:][1:I]).^2) / I)

        plot!(
            [k],
            [L2_error],
            label="",
            seriestype=:scatter,
            color=:black
        )
    end
    plot!(
        k_range,
        k_range .^ (-2),
        label="k^(-2)",
        xaxis=:log,
        yaxis=:log,
        xticks=(k_range, string.(k_range))
    )
    savefig("data/$(string(V))/optimal_diffusion_homogenized_regime_I_$(I)_d_opt_l2_error.png")

    # Plot spectral gaps convergence towards Λ_hom⋆ = 4π^2 / Z
    plot()
    for k in k_range
        a_rounded = round(a, sigdigits=2)
        b_rounded = round(b, sigdigits=2)
        dir_string = "data/" * string(V) * "/I_$(I)_a_$(a_rounded)_b_$(b_rounded)/optimal_diffusion_homogenized_regime_k_$(k)/"
        d_opt_gap = readdlm(dir_string * "d_opt_gap.txt")
        plot!(
            [k],
            d_opt_gap[:],
            label="",
            seriestype=:scatter,
            color=:black
        )
    end
    hline!(
        [4 * π^2 / Z],
        label="Lambda_hom^*",
        legend=:right
    )
    savefig("data/$(string(V))/optimal_diffusion_homogenized_regime_I_$(I)_d_opt_gap.png")

    # Plot minimum of diffusion coefficient
    plot()
    for k in k_range
        a_rounded = round(a, sigdigits=2)
        b_rounded = round(b, sigdigits=2)
        dir_string = "data/" * string(V) * "/I_$(I)_a_$(a_rounded)_b_$(b_rounded)/optimal_diffusion_homogenized_regime_k_$(k)/"
        d_opt_min = readdlm(dir_string * "d_opt_min.txt")
        plot!(
            [k],
            [d_opt_min[:]],
            label="",
            seriestype=:scatter,
            color=:black
        )
    end

    savefig("data/$(string(V))/optimal_diffusion_homogenized_regime_I_$(I)_d_opt_min.png")
end

V = sin_two_wells
k_range = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
I = 200
a = 0.0
b = Inf
optimal_diffusion_homogenized_regime(V, I, k_range, a, b)