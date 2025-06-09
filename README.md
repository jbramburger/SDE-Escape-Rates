# **Bounding Escape Rates and Approximating Quasi-Stationary Distributions of Brownian Dynamics**

This repository contains MATLAB scripts to reproduce the data and figures from [Bounding escape rates and approximating quasi-stationary distributions of Brownian dynamics](https://arxiv.org/abs/2504.00729) by [Jason J. Bramburger](https://hybrid.concordia.ca/jbrambur/).

## **Paper Abstract**
Throughout physics Brownian dynamics are used to describe the behaviour of molecular systems. When the Brownian particle is confined to a bounded domain, a particularly important question arises around determining how long it takes the particle to encounter certain regions of the boundary from which it can escape. Termed the first passage time, it sets the natural timescale of the chemical, biological, and physical processes that are described by the stochastic differential equation. Probabilistic information about the first passage time can be studied using spectral properties of the deterministic generator of the stochastic process. In this work we introduce a framework for bounding the leading eigenvalue of the generator which determines the exponential rate of escape of the particle from the domain. The method employs sum-of-squares programming to produce nearly sharp upper and lower bounds on the leading eigenvalue, while also giving good approximations of the associated leading eigenfunction, the quasi-stationary distribution of the process. To demonstrate utility, the method is applied to prototypical low-dimensional problems from the literature.

## **Required MATLAB Packages**
All scripts require YALMIP and MOSEK to run. Both packages can be download for free at 
- YALMIP: https://yalmip.github.io/download/
- MOSEK: https://www.mosek.com/downloads/

## **Repository Contents**
This repository contains MATLAB script to reproduce the results for the examples in Section 4 of the paper. Precisely, the scripts perform the following tasks:
- `double_well.m` applies the method to the escape from a double well potential landscape and reproduces all results presented in Section 4.1
- `ball_escape_2d.m` applies the method to the escape from the unit ball in 2D and reproduces results presented in Section 4.2.1
- `ball_escape_3d.m` applies the method to the escape from the unit ball in 3D and reproduces results presented in Section 4.2.1
- `ball_escape_2_holes.m` applies the method to the escape from the unit ball in 2D with 2 symmetrically placed small exits and reproduces results presented in Section 4.2.2
- `ball_escape_1_hole.m` applies the method to the escape from the unit ball in 2D with 1 small exit and reproduces results presented in Section 4.2.2
