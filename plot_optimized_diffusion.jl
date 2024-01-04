"""Plot the optimal diffusion obtained from the optimal_diffusion_computation.jl script.
"""

using Plots, DelimitedFiles
include("potentials.jl")

I = 1000
a = 0.0
b = Inf
p = 2
a_rounded = round(a, sigdigits = 2)
b_rounded = round(b, sigdigits = 2)
V = sin_two_wells
dir_string = "data/" * string(V) * "/I_$(I)_a_$(a_rounded)_b_$(b_rounded)/"
d_opt = readdlm(dir_string * "d_opt.txt")
d_homog = map(x -> exp(V(x)), XX)

function lp_constraint(d, I, mu_arr, p)
    return (dot(d .^ p, mu_arr .^ p) / I)^(1 / p)
end
XX = [i / I for i = 0:(I-1)]
mu_arr = map(x -> exp(-V(x)), XX)
d_constant = ones(I) / lp_constraint(ones(I), I, mu_arr, p)

# potential
plot(XX, V.(XX), label = "potential energy function", color = :green, linewidth = 2)
savefig(dir_string * "potential_energy_function.png")

plot(XX, d_opt, label = "optimal diffusion", color = :red, linewidth = 2)
plot!(
    XX,
    d_homog,
    label = "homogenized diffusion",
    color = :blue,
    linestyle = :dash,
    linewidth = 2,
)
plot!(
    XX,
    d_constant,
    label = "constant diffusion",
    color = :black,
    linestyle = :dashdot,
    linewidth = 2,
)
savefig(dir_string * "d_opt.png");
