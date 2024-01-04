"""One-dimensional optimization procedure to compute numerically the optimal diffusion for the overdamped Langevin process on the torus.

This script uses IPOPT for the optimization procedure.
"""

using LinearAlgebra, JuMP, Ipopt, DelimitedFiles

# Variables optimized in the code are x[i] =  exp(-\beta V(q_i))D(q_i)

# Generalized eigenvalue problem: Au=\lambda Mu

"""
    optim_algo(
    V,
    I,
    pi_arr::AbstractVector{T},
    Z;
    p=2.0,
    a=0.0,
    b=Inf,
    tol=10^(-2),
    max_it=200,
    save=true,
    rewrite_save=true
) where {T}

Optimization algorithm to find optimal diffusion

# Arguments
- `n::Integer`: the number of elements to compute.
- `V::Function`: potential energy function defined on the torus
- `I::Int`: number of hat functions in P1 Finite Elements basis
- `pi_arr::AbstractVector{T}`: discrete approximation of the Boltzmann-Gibbs measure
- `Z::Real`: partition function, equal to int exp(-V)
- `p::Real`: L^p constraint. Defaults to 2..
- `a::Real`: lower bound for the variable constraint. Defaults to 0..
- `b::Real`: upper bound for the variable constraint. Defaults to Inf.
- `tol::Real`: overall NLP error for the IPOPT algorithm. Defaults to 10^(-2).
- `max_it`: maximum number of iterations for the IPOPT algorithm. Defaults to 200.
- `save::Bool`: if saving first, second, third and fourth eigenvalues and the constraint and gradient norm values during the optimization procedure. Defaults to true.
- `rewrite_save::Bool`: if forcing save by rewriting previously saved data. Defaults to true.
"""
function optim_algo(
    V,
    I,
    mu_arr::AbstractVector{T};
    p = 2.0,
    a = 0.0,
    b = Inf,
    tol = 10^(-2),
    max_it = 200,
    save = true,
    rewrite_save = true,
) where {T}

    """
        lp_constraint(
            x, I, p
        )

        Compute the L^p constraint for x
    """
    function lp_constraint(x, I, p)
        return (sum(x .^ p) / I)^(1 / p)
    end

    """
        construct_M(I, mu_arr)

        Construct the matrix appearing on the right hand side of the generalized eigenvalue problem
    """
    function construct_M(I, mu_arr::AbstractVector{T}) where {T}
        M = zeros(T, I, I)
        for i = 1:I-1
            M[i, i] += mu_arr[i] / (3 * I)
            M[i+1, i+1] += mu_arr[i] / (3 * I)
            M[i, i+1] += mu_arr[i] / (6 * I)
            M[i+1, i] += mu_arr[i] / (6 * I)
        end
        M[I, I] += mu_arr[I] / (3 * I)
        M[1, 1] += mu_arr[I] / (3 * I)
        M[I, 1] += mu_arr[I] / (6 * I)
        M[1, I] += mu_arr[I] / (6 * I)
        return M
    end

    """
        f(x::T...) where {T<:Real}

        Objective function to be maximized. Outputs the spectral gap.
    """
    function f(x::T...) where {T<:Real}

        # Construct the matrix depending on D
        A = zeros(I, I)
        for i = 1:I-1
            A[i, i] += x[i] * I
            A[i+1, i+1] += x[i] * I
            A[i, i+1] -= x[i] * I
            A[i+1, i] -= x[i] * I
        end
        A[I, I] += x[I] * I
        A[1, 1] += x[I] * I
        A[I, 1] -= x[I] * I
        A[1, I] -= x[I] * I

        # Obtain eigenvalues
        eig = eigen(A, M)
        vals = eig.values
        idx = sortperm(vals)
        vals = vals[idx]

        return vals[2]
    end

    """
        ∇f(g::AbstractVector{T}, x::T...) where {T<:Real}

        Gradient of the objective function. g is filled in place.
    """
    function ∇f(g::AbstractVector{T}, x::T...) where {T<:Real}

        # Construct the matrix depending on D
        A = zeros(I, I)
        for i = 1:I-1
            A[i, i] += x[i] * I
            A[i+1, i+1] += x[i] * I
            A[i, i+1] -= x[i] * I
            A[i+1, i] -= x[i] * I
        end
        A[I, I] += x[I] * I
        A[1, 1] += x[I] * I
        A[I, 1] -= x[I] * I
        A[1, I] -= x[I] * I

        # Obtain eigenvalues and eigenvectors
        eig = eigen(A, M)
        vals = eig.values
        vecs = eig.vectors
        idx = sortperm(vals)
        vals = vals[idx]
        val = vals[2]
        vecs = vecs[:, idx]

        # for ergonomic reasons for gradient computation
        vp = zeros(T, I + 1)
        vp[2:I+1] .= vecs[:, 2]
        vp[1] = vp[I+1]

        # matrix corresponding to the tri-diagonal matrix A
        norm_factor = vp[1:I]' * M * vp[1:I]

        for i = 1:I
            g[i] = (vp[i+1]-vp[i])^2 * I / norm_factor
        end

        if save
            push!(first_eig, vals[1])
            push!(second_eig, val)
            push!(third_eig, vals[3])
            push!(fourth_eig, vals[4])
            #push!(norm_jac, norm(g) / sqrt(I))
            push!(constraints, lp_constraint(x, I, p))
        end
        return
    end

    function h(x::T...) where {T<:Real}
        return sum(x .^ p) / I
    end

    function ∇h(g::AbstractVector{T}, x::T...) where {T<:Real}
        for i = 1:I
            g[i] = p * x[i]^(p - 1) / I
        end
        return
    end


    if save
        a_rounded = round(a, sigdigits = 2)
        b_rounded = round(b, sigdigits = 2)
        dir_string = "data/" * string(V) * "/I_$(I)_a_$(a_rounded)_b_$(b_rounded)/"

        if !rewrite_save
            if isfile(dir_string * "first_eigenvalue.txt")
                return
            end
        end

        first_eig = []
        second_eig = []
        third_eig = []
        fourth_eig = []
        norm_jac = []
        constraints = []
    end

    # Construct M once and for all
    M = construct_M(I, mu_arr)

    # First guess for the optimization procedure, D=1/mu (homogenized diffusion)
    x_init = ones(I)

    # Create optimiation model using the IPOPT optimizer
    model = Model(Ipopt.Optimizer)
    # Set a maximum number of iterations to 2000
    set_optimizer_attribute(model, "max_iter", max_it)
    # Set tolerance for the overall NLP error
    set_optimizer_attribute(model, "tol", tol)

    # Register user defined functions
    register(model, :my_f, I, f, ∇f; autodiff = false)
    register(model, :my_h, I, h, ∇h; autodiff = false)

    # Variable to be optimized, equal to D*mu
    @variable(model, x[1:I])

    # Set first guess
    set_start_value.(x, x_init)

    # Set objective value
    @NLobjective(model, Max, my_f(x...))

    # Set constraints
    @NLconstraint(model, my_h(x...) <= 1) # L^p constraint
    @NLconstraint(model, [i = 1:I], a <= x[i] <= b) # bounds constraints

    print(model)
    optimize!(model)

    # If failure
    term_status = termination_status(model)
    if term_status != LOCALLY_SOLVED && term_status != ALMOST_LOCALLY_SOLVED
        println(term_status)
        println("Not saving")
        return
    end

    if term_status == ALMOST_LOCALLY_SOLVED
        println(term_status)
        println("Beware, this was only almost locally solved")
    end

    x_opt = value.(x)
    XX = [i / I for i = 0:(I-1)]
    inv_mu_arr = map(x -> exp(V(x)), XX)
    d_opt = x_opt .* inv_mu_arr
    gap_opt = objective_value(model)
    # If success, saves
    if save
        min_d = minimum(d_opt)

        mkpath(dir_string)

        writedlm(dir_string * "first_eigenvalue.txt", first_eig)
        writedlm(dir_string * "second_eigenvalue.txt", second_eig)
        writedlm(dir_string * "third_eigenvalue.txt", third_eig)
        writedlm(dir_string * "fourth_eigenvalue.txt", fourth_eig)
        writedlm(dir_string * "d_opt.txt", d_opt)
        writedlm(dir_string * "d_opt_gap.txt", [gap_opt])
        writedlm(dir_string * "d_opt_min.txt", [min_d])
        writedlm(dir_string * "constraint.txt", constraints)
        writedlm(dir_string * "gradient_norm.txt", norm_jac)
    end

    return d_opt, gap_opt
end
