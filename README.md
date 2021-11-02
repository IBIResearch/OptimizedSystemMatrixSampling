# Optimized System Matrix Sampling Example

This folder contains example code for the generation of optimized sampling patterns for system matrices in magnetic particle imaging (MPI).
The method uses the sparsity structure of a previously measured system matrix and searches for an optimized sampling pattern based on the Bayesian Fisher information matrix. The resulting patterns can then be used for the compressed sensing based recovery of other system matrices with similar structure.

The method is described in the associated publication

M. Grosser, T. Knopp, Optimized sampling patterns for the sparse recovery of system matrices in Magnetic Particle Imaging, *International Journal on Magnetic Particle Imaging* <!--, 2021  [*arXiv:2101.12624*](https://arxiv.org/abs/2006.05741). -->


## Installation

In order to use this code one first has to download [Julia](https://julialang.org/) (version 1.6 or later), clone this repository and navigate to the folder in the command line. The example script automatically activates the environment and installs all necessary packages.

## Execution
After installation the example code can be executed by running `julia` and entering
```julia
include("example.jl")
```
This will first download all data and then compute an optimized sampling pattern based on a 2d system matrix of liquid Perimag particles. Subsequently, matrix recovery of a 2d system matrix of immobilized Perimag particles is performed.
Parameters of the pattern optimization and matrix recovery are documented in the Julia script and can be changed. After computation of the optimized sampling patterns, the script will open a plotting window and will show the obtained patterns together with Poisson disk patterns at the same undersampling factors. Similarly, the script will open a plotting window and show some frequency components of the recovered system matrices after matrix recovery is performed.

## Open MPI Data

The measurement data associated to this project is about 33 MB large and will be downloaded and stored automatically, when the code is executed for the first time.
It is published under a [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/legalcode) license and can be found here:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5638874.svg)](https://doi.org/10.5281/zenodo.5638874)
