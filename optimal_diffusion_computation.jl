include("main.jl")
include("potentials.jl")

I = 500
p = 2
a = 0.0
b = Inf
V = sin_two_wells
tol = 10^(-1)

function mu(x)
    return exp(-V(x))
end

XX = [i / I for i in range(0, I)]
mu_arr = map(x -> exp(-V(x)), XX)
Z = sum(mu_arr) / I
pi_arr = mu_arr / Z

optim_algo(
    V, I, pi_arr, Z,
    p=p,
    a=a, b=b,
    tol=tol
);