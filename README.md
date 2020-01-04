# Spiral-Waves-Boundary-Sinks-and-Spectra
Matlab code to reproduce computations in [Dodson and Sandstede 2019].

# Instructions for Use

This repository provides Matlab code to solve for spiral wave patterns, spectra, and additional related tasks from [S Dodson and B Sandstede. Determining the source of period-doubling instabilities in spiral waves. SIAM Journal on Applied Dynamical Systems 18 (2019) 2202-2226](http://dx.doi.org/10.1137/19M1264813). Details about the methods can be found in the paper, here we provide instructions on how to run the codes.

In general, the codes are organized into folders by task, and the contents of each folder are described in the relevant sections. The folder `Util` contains utility functions that one or more scripts may call. The scripts work for the Karma and R&ouml;ssler models by calling the model specific functions, which lie in the `Util` folder, and much of this code can be adapted for other reaction-diffusion systems by modifying these model specific functions.

## Data Files

Contains files for R&ouml;ssler and Karma model asymptotic wave trains, spirals on bounded disks, and boundary sinks which can be used as initial conditions for solving or continuing systems. Files are labeled with model and key system or numerical parameters.

## Solve Patterns

All patterns are set up as root finding problems using the built-in Matlab solver `fsolve`. The Matlab scripts used to solve each of the patterns is listed below.

**Asymptotic Wave Trains:**  `solve_1D_wave_train`

**Spiral Wave on Disks:**  `solve_spiral_pattern`
Works for spiral on disks with Neumann boundary conditions or non-reflecting boundary conditions, depending on functions specified inthe variable `file_names.solver_fcn` and the phase condition given by
`file_names.set_up_phase_condition`.

**Boundary Sink:**  `solve_boundary_sink`
Methods based on those in [Lloyd and Scheel, 2017](https://epubs.siam.org/doi/10.1137/16M1073212) and [Goh and Scheel, 2018](https://londmathsoc.onlinelibrary.wiley.com/doi/full/10.1112/jlms.12122).

## Continue Patterns

The continuation code uses a mixture of secant and simple continuation methods. 

**Asymptotic Wave Trains:** `continue_wave_train`
Secant continuation.

**Spiral Waves on Disks:** `continue_spiral`
Secant continuation that works for spiral waves on bounded disks with Neumann or non-reflecting boundary conditions based on solver files. Continue from a Neumann boundary condition to non-reflecting using the non-reflecting options (`file_names.problem ='Rossler_2D_spiral_nonreflecting'` and
`file_names.set_up_phase_condition ='spiral_non_reflecting_phase_condition'`) with continuation parameter &kappa; (`contPar.Name = 'kappa'`).

**Point Eigenfunction:** `continue_non_reflecting_bc_eigenfunction`
Simple continuation that starts with spiral and eigenfucntion with Neumann boundary conditions and continues to non-reflecting boundary conditions using the two-step procees of (1) solving for spiral wave and (2) computing new eigenfunction. The two-step process is necessary because the linearization in the eigenvalue problem requires a spiral solution with the same parameter values.

## Spectra

Functions to compute the spectra of spirals on bounded disks and boundary sink.

**Essential Spectrum:** `compute_essential_spectrum`
Simple continuation to compute the essential spectrum. Methods follow those outlined in [Rademacher et al 2007](https://www.sciencedirect.com/science/article/pii/S0167278907000863?via%3Dihub).

**Absolute Spectrum:** `continue_absolute_spectrum`
Simple continuation along absolute spectrum curve. An initial starting point is provided for the Karma and Rossler systems in the `data_files/` folder. Methods follow those outlined in [](Rademacher et al 2007). The function defined by `file_names.spatial_evals_fcn = 'spatial_evals_fcn_karma'` computes the spatial eigenvalues &nu; for the given value of the temporal eigenvalue &lambda;. By operating independently from the continuation this function provides a second check of the spatial eigenvalue distribution.

**Spiral Point Spectrum:** `compute_point_spectra`
Uses the sparse eigenvalue solver `eigs` to compute eigenvalues of spiral on a bounded disk near the user defined seed points (`seed_real` and `seed_imag`). Works for Neumann or non-reflecting boundary conditions based on defined jacobian fuction. Warning: do not use `eigs` if the system is highly non-normal, especially if one of the diffusion coefficients is 0.

**Boundary Sink Point Spectrum:** `compute_point_spectra_boundary_sink`
Uses the sparse eigenvalue solver `eigs` to compute eigenvalues of spiral on a bounded disk near the user defined seed points (`seed_real` and `seed_imag`).

