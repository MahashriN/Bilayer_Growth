# Bilayer Growth

This repository contains MATLAB code used to generate numerical results for the paper

*“Heterogeneous growth enables pattern formation in coupled bilayer reaction–diffusion systems.”*

The simulations study pattern formation in two coupled layers on a one-dimensional spatial domain, incorporating diffusion, nonlinear reaction kinetics, and growth-induced effects.

The numerical implementation is based on a finite element discretisation in one spatial dimension.

### How to run the code

1. Run the main simulation

   Execute
   `Main2_1D2L_Diff_Iso.m`

   This script:

   - defines all model parameters and initial conditions,
   - assembles the finite element matrices,
   - evaluates reaction kinetics and inter-layer coupling,
   - advances the coupled bilayer growing system in time,
   - and produces the numerical output used in the paper.

### Main functions used

The main script calls the following core routines:

 - `AssembleGlobalMatrices1D.m` - Assembles the global mass and stiffness matrices for the one-dimensional finite element discretisation.
 - `ReactKineInt1D.m` - Computes the nonlinear reaction kinetics in one dimension.
 - `SA2_1D2L_Diff_Iso.m` - Sets up the frozen linear stability problem for the coupled bilayer recton-diffusion system and computes the corresponding dynamic dispersion relation.

These routines rely on:

 - `basis_linear_1D.m` - Linear basis functions on the one-dimensional reference element.
 - `RefEdgeQuad.m` - Reference quadrature definitions used for numerical integration.

### Visualisation utilities

 - `cmap_colorbar.m` - Custom colormap and colorbar utility for visualising spatiotemporal patterns.
 - `brownyellow.mat` - Stored colormap data used by the plotting routines.

### Notes
- All parameters controlling diffusion, coupling, kinetics, and growth effects are defined inside the main script.
- The code is written for clarity and reproducibility rather than computational optimisation.
- This repository is intended to accompany the above paper and reproduce its numerical results.
