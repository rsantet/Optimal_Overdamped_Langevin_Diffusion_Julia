using DelimitedFiles, LinearAlgebra
include("potentials.jl")

"""
    d_interp(d, I, x)

Outputs the value of the diffusion coefficient d at point x using a piecewise linear approximation

# Arguments
- `d::Vector`: the diffusion coefficient values for the point in the mesh
- `I::Int`: number of point in the mesh (length of vector d)
- `x0::Real`: position where to compute the diffusion coefficient
"""
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


function get_one_transition(d, I, Δt, x0, V)
    i = 0
    x = x0
    sqrt_2_Δt = sqrt(2 * Δt)
    MH_counter = 0
    while x0 - 1 <= x <= x0 + 1
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
            MH_counter += 1
        else
            x = proposal
        end
        i += 1
    end
    return i * Δt, MH_counter / i
end


function transition_times(d, I, Δt, N, x0, V)

    times = []
    MH = []
    for i = 0:N-1
        if i % 100 == 0
            println("Iteration $(i)/$(N)")
            flush(stdout)
        end

        T, MH_ratio = get_one_transition(d, I, Δt, x0, V)
        push!(times, T)
        push!(MH, MH_ratio)

    end

    return times, MH
end



V = sin_two_wells
# global min of potential is 0.3654418277735119 #
# local min of potential is 0.8971356821573657 #
x0 = 0.3654418277735119
N = 10000

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
d_homog = 1.0 ./ mu_arr

# Constant diffusion
d_constant = ones(I) / lp_constraint(ones(I), I, mu_arr, p)

Δt_range = [0.00005, 0.0001, 0.0005, 0.001, 0.005, 0.01]
for (idx, Δt) in enumerate(Δt_range)
    for (d, d_string) in
        zip([d_opt, d_homog, d_constant], ["d_opt", "d_homog", "d_constant"])
        file_string =
            ("data/" * string(V) * "/transition_times_Δt_$(Δt)_N_$(N)_$(d_string)")
        if !isfile(file_string * ".txt")
            println("Running for Δt = $(Δt) and $(d_string)")
            flush(stdout)
            times, MH_ratios = transition_times(d, I, Δt, N, x0, V)
            writedlm(file_string * ".txt", times)
            writedlm(file_string * "_MH_ratios.txt", MH_ratios)
        end
    end
end

# For various lower bounds
a_range = [0.05 * i for i = 0:20]
Δt_range = [0.00005, 0.0001, 0.0005, 0.001, 0.005, 0.01]
for a in a_range
    local a_rounded = round(a, sigdigits = 2)
    dir_string = "data/" * string(V) * "/I_$(I)_a_$(a_rounded)_b_$(b_rounded)/"
    d_opt = readdlm(dir_string * "d_opt.txt")

    for Δt in Δt_range
        file_string =
            ("data/" * string(V) * "/transition_times_Δt_$(Δt)_N_$(N)_d_opt_a_$(a_rounded)")
        if !isfile(file_string * ".txt")
            println("Running for Δt = $(Δt) and a = $(a_rounded)")
            flush(stdout)
            times, MH_ratios = transition_times(d_opt, I, Δt, N, x0, V)
            writedlm(file_string * ".txt", times)
            writedlm(file_string * "_MH_ratios.txt", MH_ratios)
        end
    end
end