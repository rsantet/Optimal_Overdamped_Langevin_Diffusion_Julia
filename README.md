# Optimal Overdamped Langevin diffusion

This repository contains Julia code to optimize the diffusion of overdamped Langevin dynamics. The method is described in the article "Optimizing the diffusion of overdamped Langevin dynamics" by T. Lelièvre, G. Pavliotis, G. Robin, R. Santet and G. Stoltz.

## Method

Overdamped Langevin dynamics are used in many MCMC sampling algorithms to produce new states. An important parameter of such dynamics is the diffusion function which has fundamental impact on the convergence rate of Langevin dynamics. The method developed in this repository aims to compute the optimal diffusion function to maximize the convergence speed of overdamped Langevin dynamics towards the target measure. The resulting optimal diffusion can be plugged into MCMC algorithms (RWMH, MALA, etc.) to accelerate their convergence towards equilibrium. The code available here concerns only *one-dimensional* problems on the *torus*, but can be adapted to other low-dimensional settings.

## Files

The repository contains:
- main.jl: the optimization algorithm
- potentials.jl: examples of one-dimensional potentials (minus log of target measure)
- optimal_diffusion_computation.jl: script to perform the optimization algorithm on one-dimensional potentials
- plot_optimized_diffusion.jl: plot the optimized diffusion obtained from the main.py script
- RWMH.jl: a Random Walk Metropolis-Hastings algorithm implementation where the optimal diffusion computed with optimal_diffusion_computation.jl can be used as a parameter
- plot_RWMH.jl: plot (unperiodized) trajectories and sampling results of the RWMH algorithm using either the optimal, homogenized or constant diffusion coefficient
- convergence_RWMH.jl: compute the TV and $L^2$ distance between the empirical distribution of the MC and the equilibrium distribution for a 1D double well potential.
- plot_convergence_RWMH.jl: plot the results of the previous script
- msd.jl: compute and plot the Mean Square Displacement of the MCs induced by various diffusion coefficients for a 1D double well potential.
- transition_times.jl: compute the duration time to transition between two metastable states of a 1D double well potential for various diffusion coefficients (optimal, homogenized, constant)
- plot_transition_times.jl: plot the results obtained from the previous script
- optimal_diffusion_homogenized_regime.jl: visualization of the behaviour of the optimal diffusion in the asymptotic regime of periodic homogenization. In particular, it shows that (1) the optimal diffusion is $1/k$-periodic, (2) the sequence of rescaled optimal diffusion converges in $L^2$ towards $D_{\mathrm{hom}}^\star$ at rate $\mathrm{O}(k^{-2})$, (3) the sequence of spectral gaps converges to $\Lambda_{\mathrm{hom}}^\star=4\pi^2 / Z$ where $Z=\int_{\mathbb{T}}\mathrm{e}^{-V}$, (4) the minimum of the optimal diffusion quickly goes away from zero
- varying_lower_bound.jl: varying_lower_bound.py: compute the optimal diffusion coefficient for various lower bounds $a\in[0,1]$
- rejection_probability_spatial_analysis.jl: plot the rejection probability as a function of space for MCs induced by various diffusion coefficients for a 1D double well potential.

## TODO
- Add a comparison with smooth-min approach with analytic gradient formula