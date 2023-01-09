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
    RWMH(d, I, Δt, N_it, x0, V, β, Gs=[])

Performs a Random Walk Metropolis-Hastings algorithm using the diffusion coefficient d

# Arguments
- `d::Vector`: the diffusion coefficient values for the point in the mesh
- `I::Int`: number of point in the mesh (length of vector d)
- `Δt::Real`: time step for the simulation
- `N_it::Int`: number of time steps for the simulation
- `x0::Real`: initial position
- `V::Function`: potential energy function defined on the torus
- `β::Real`: inverse temperature
- `Gs::Vector`: Optional. The gaussian increments used for the simulation. Must be of length N_it. Default to randn(N_it).
"""
function RWMH(d, I, Δt, N_it, x0, V, β, Gs=[])
    trajectory = []
    x = x0
    sqrt_2_dt = sqrt(2 * Δt / β)
    
    if length(Gs) != N_it
        Gs = randn(N_it)
    end

    MH_counter = 0
    for i in 1:N_it
        
        if i%10000 == 0
            println("$(i)/$(N_it)")
            flush(stdout)
        end

        G = Gs[i]
        d_x = d_interp(d, I, x)
        proposal = x + sqrt_2_dt * sqrt(d_x) * G
        d_proposal = d_interp(d, I, proposal)
        V_x = V(x)
        V_proposal = V(proposal)
        sqrt_d_ratio = sqrt(d_x / d_proposal)
        G_proposal = sqrt_d_ratio * G
        alpha = log(sqrt_d_ratio) - β*(V_proposal - V_x) - (G_proposal^2 - G^2) / 2

        if log(rand()) > alpha
            MH_counter += 1
        else
            x = proposal
        end
        push!(trajectory, x)
    end

    MH_rejection_probability = MH_counter / N_it
    return trajectory, MH_rejection_probability
end