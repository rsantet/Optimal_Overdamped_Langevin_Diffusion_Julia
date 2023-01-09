"""Plot the optimal diffusion obtained from the optimal_diffusion_computation.jl script.
"""

using Plots, DelimitedFiles
include("potentials.jl")

I = 500
a = 0.0
b = Inf
a_rounded = round(a, sigdigits=2)
b_rounded = round(b, sigdigits=2)
XX = [i / I for i in 0:(I-1)]
V = sin_two_wells
dir_string = "data/" * string(V) * "/I_$(I)_a_$(a_rounded)_b_$(b_rounded)/"
d_opt = readdlm(dir_string * "d_opt.txt")
d_homog = map(x -> exp(V(x)), XX)

plot(XX, d_opt, label="Optimal diffusion")
plot!(XX, d_homog, label="Homogenized diffusion")
savefig(dir_string * "d_opt.png");