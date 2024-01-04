include("main.jl")
include("potentials.jl")

I = 1000
p = 2
a = 0.0
b = Inf
V = sin_two_wells
tol = 10^(-12)

XX = [i / I for i = 0:I-1]
mu_arr = map(x -> exp(-V(x)), XX)
Z = sum(mu_arr) / I
pi_arr = mu_arr / Z

optim_algo(V, I, mu_arr, p = p, a = a, b = b, tol = tol);
