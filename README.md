# Bilayer Growth

This repository contains MATLAB code used to simulate pattern formation in coupled reaction–diffusion systems on growing domains.

The code investigates pattern formation in two coupled reaction–diffusion layers on a one-dimensional growing domain. The model incorporates diffusion, nonlinear reaction kinetics, inter-layer coupling, and growth-induced dilution and geometric effects.

The numerical implementation is based on a one-dimensional finite element discretisation coupled with an IMEX time-stepping scheme.

## How to run the code

### Main simulation

Run

`main_simulation.m`

This script:

* defines all model parameters and initial conditions,
* assembles the finite element matrices,
* evaluates reaction kinetics and inter-layer coupling,
* advances the coupled bilayer system in time,
* and generates the numerical solutions and figures.

### Frozen-time stability analysis

Run

`frozen_time_stability_analysis.m`

This script performs the frozen-time linear stability analysis of the coupled bilayer system and computes the corresponding dynamic dispersion relations used to investigate pattern-forming instabilities.

## Main functions

The main simulation relies on the following routines:

* `AssembleGlobalMatrices1D.m` – assembles the global mass and stiffness matrices.
* `ReactKineInt1D.m` – evaluates the nonlinear reaction kinetics.

These routines make use of:

* `basis_linear_1D.m` – linear basis functions on the reference element.
* `RefEdgeQuad.m` – quadrature rules used for numerical integration.

## Visualisation utilities

* `cmap_colorbar.m` – custom colormap and colorbar utility.
* `brownyellow.mat` – stored colormap data used by the plotting routines.

## Reproducibility

All parameters controlling diffusion, reaction kinetics, inter-layer coupling, and domain growth are specified directly within the main scripts.

The code is intended to reproduce the numerical experiments reported in the accompanying manuscript and is written with clarity and reproducibility in mind.
