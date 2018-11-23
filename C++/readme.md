# Table of contents

## 1. [New_simulation](#it1)
## 2. [Maths](#it2)
## 3. [EDWA_sim](#it3)
## 4. [spot_size_calculator](#it4)

## New_simulation <a name = "it1"></a>
This is a finite difference time domain method solver for Maxwell's Equations, featuring second order absorbing boundary
in 1D and 2D grids, also the ability to easily add Perfect electrical conductors (PEC) and dielectric structures and 
be able to choose from a variety of sources.

<b>Note:</b> the text files beginning with plot_script_  are GNUPLOT scripts used to quickly visualize the data
produced by the simulation

## Maths <a name = "it2"></a>
Contains a set of useful headers defining a Matrix class with error handling and more robust algorithms. It also has defined some parallel operations for matrices 
with a total amount of elements greater than 500k. It contains the following operations:
* Matrix slicing
* Matrix transpose
* Matrix invert
* Matrix determinant (via Cholesky or LU)
* Matrix Eigenvalues
* LU (Lower upper) factorization
* Cholesky Factorization
* QR factorization (Gramm Schmidt orthonormalization process)
Library dependencies: pthread

## EDWA_sim (in progress) <a name = "it3"></a>
Four level Erbium, Two level Ytterbium simulation which takes into account up-conversion, cross-relaxation between Ytterbium and Erbium. Currently in progress. GUI made using [dear imgui](https://github.com/ocornut/imgui) with [SFML bindings](https://github.com/eliasdaler/imgui-sfml). Added support for
concurrent execution of multiple simulations with progress bar for each thread of execution. Requires [gnuplot](http://www.gnuplot.info/) installed
Library dependencies: pthread, OpenGL, sfml-graphics, sfml-window, sfml-system, boost_system, boost_filesystem


## Spot_size_calculator <a name = "it4"></a>
This program takes the pixels/length and an image of the laser spot with a dark background (ideally in a dark room) as input and bulk processes all laser spot image 
files provided, adding up the number of pixels with an intensity value >= ~13.5% of max and then multiplying that value by the area of each pixel (assuming square pixels) 
to obtain the spot size. The units will depend on the units given when entering the pixels/length value (e.g. if pixels/mm given then area will be mm^2)
Library dependencies: 

