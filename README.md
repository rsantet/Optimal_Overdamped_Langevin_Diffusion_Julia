# Optimal Overdamped Langevin diffusion

This repository contains Julia code to optimize the diffusion of overdamped Langevin dynamics. The method is described in the article "Optimizing the diffusion of overdamped Langevin dynamics" by T. Lelièvre, G. Pavliotis, G. Robin, R. Santet and G. Stoltz.

## Method

Overdamped Langevin dynamics are used in many MCMC sampling algorithms to produce new states. An important parameter of such dynamics is the diffusion function which has fundamental impact on the convergence rate of Langevin dynamics. The method developed in this repository aims to compute the optimal diffusion function to maximize the convergence speed of overdamped Langevin dynamics towards the target measure. The resulting optimal diffusion can be plugged into MCMC algorithms (RWMH, MALA, etc.) to accelerate their convergence towards equilibrium. The code available here concerns only *one-dimensional* problems on the *torus*, but can be adapted to other low-dimensional settings.

## Files

The repository contains:
- main.jl: the optimization algorithm
- potentials.jl: examples of one-dimensional potentials (minus log of target measure)
- RWMH.jl: a Random Walk Metropolis-Hastings algorithm implementation where the optimal diffusion computed with main.jl can be used as a parameter
- plot_optimized_diffusion.jl: plot the optimized diffusion obtain from the main.py script

## TODO
- Propose a visualization of the behaviour of the optimal diffusion in the asymptotic regime of periodic homogenization
- Add a comparison with smooth-min approach with analytic gradient formula