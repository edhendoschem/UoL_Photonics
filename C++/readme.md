# Table of contents
## 1. [FDTD_parallel](#FDTD_parallel)
## 2. [new_simulation](#new_simulation)

## FDTD_parallel <a name = "FDTD_parallel"></a>
This is a finite difference time domain method solver for Maxwell's Equations, to simulate propagation through a
rectangular straight waveguide. It employs multithreading to enhance speed also provides a serial execution functions
to compare speedup and correctness.
<b>Currently unfinished</b>

## new_simulation <a name = "new_simulation"></a>
This is an improved version of FDTD_parallel, featuring second order absorbing boundary in 1D and 2D grids, also
the ability to easily add Perfect electrical conductors (PEC) and dielectric structures.

<b>Note:</b> the text files plot_script, plot_script_2d, etc are GNUPLOT scripts used to quickly visualize the data
produced by the simulation

###<b>To do list:<\b>
1. Fix TFSF area leakage in corners with dispersion compensation techniques
2. Implement 3D simulation
3. Implement Parallel processing in 3D simulation
4. Implement amplification in 3D simulation
5. Implement amplification in 2D simulation
6. Implement Fourier analysis
