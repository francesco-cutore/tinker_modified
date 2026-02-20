# tinker_modified

This repository contains the source code developed for the MSc Thesis: **"X-ray Scattering-restrained molecular dynamics for the structural characterization of complex fluids and materials"** at Politecnico di Milano. 

The complete modification of the Tinker source code developed for this thesis encompasses over 1,500 lines of Fortran 95. To ensure reproducibility, ease of access, and a concise printed manuscript, the full repository has been made publicly available here.

## Overview
Molecular Dynamics (MD) force fields sometimes struggle to accurately reproduce the structural properties of complex systems, such as ionic liquids. While methods like Empirical Potential Structure Refinement (EPSR) exist, they inherently fail to account for the true time-evolution of the system.

This project introduces a method to incorporate experimental X-ray Scattering data directly into MD simulations on-the-fly. By introducing an additional harmonic biasing term into the MD potential, the simulation generates fictitious forces that guide the system into configurations consistent with the experimental structure factor $S(q)$.To optimize computational efficiency, the structure factor is evaluated via a direct calculation method rather than the traditional, computationally heavy Radial Distribution Function (RDF) approach.

## Key Features
* **Direct $S(q)$ Calculation:** Calculates the structure factor directly using real and imaginary phase components, $E(q)$ and $G(q)$, avoiding $O(N^2)$ distance calculations.
* **On-the-Fly Restraint Forces:** Derives analytical gradients to distribute harmonic restraint forces to each atom at every integration step.
* **Discrete Wavevector Generation:** Generates a discrete grid of wavevectors perfectly matching Periodic Boundary Conditions (PBC) to prevent boundary artifacts.
* **High-q Decimation:** Optimizes performance by applying a power-law decimation algorithm to reduce the density of sampled vectors in the computationally expensive high-q region.
* **Precomputed Form Factors:** Computes X-ray atomic form factors and normalization denominators during initialization to minimize per-step overhead.
* **OpenMP Parallelization:** Accelerates heavy array operations for the force calculations to maintain simulation efficiency.

## Core Modules & Subroutines
This implementation modifies the open-source **TINKER** molecular modeling package. The core additions include:
* `radialask.f`: The global controller that initializes the wavevector grid, maps atom types and mole fractions, and loads experimental target data.
* `generate_lattice.f`: Generates the discrete grid of $q$-vectors based on the exact simulation box dimensions.
* `decimation.f`: Reduces the number of high-q wavevectors to speed up the simulation.
* `precompute_form_factors.f`: Calculates and caches atomic form factors and normalization constants.
* `structfactor_forces.f`: The main routine called at each MD time step. It calculates the instantaneous structure factor, evaluates the harmonic penalty energy, and distributes the restraint forces back to the atoms.

## Validation Systems
The methodology and algorithms were benchmarked and validated against simple fluids:
1. **Liquid Argon (94.4 K):** Used to test the pure effect of the restraint potential on Lennard-Jones interactions, demonstrating the ability to drive structural evolution. 
2. **Liquid Water (TIP3P, 298 K):** Validated the system-size scaling and minimum wavevector limits, comparing the direct calculation method against TRAVIS RDF-derived structure factors.

## Author & Acknowledgments
* **Author:** Francesco Cutore 
* **Advisors:** Prof. Guido Raos, Dott. Alessandro Mariani [cite: 61, 62]
* **Institution:** School of Industrial and Information Engineering, Politecnico di Milano [cite: 53]
* **Academic Year:** 2025-26 