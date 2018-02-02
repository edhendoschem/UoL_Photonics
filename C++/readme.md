# Table of contents

## 1. [new_simulation](#it1)
## 2. [useful_headers](#it2)

## new_simulation <a name = "it1"></a>
This is a finite difference time domain method solver for Maxwell's Equations, featuring second order absorbing boundary
in 1D and 2D grids, also the ability to easily add Perfect electrical conductors (PEC) and dielectric structures and 
be able to choose from a variety of sources.

<b>Note:</b> the text files beginning with plot_script_  are GNUPLOT scripts used to quickly visualize the data
produced by the simulation

## useful_headers <a name = "it2"></a>
Contains:

- aux_maths: Contains Monte Carlo integration, bisection method, Runge Kutta method, useful constants 
- Flat_vec: A vector with defined access operator in order to be able to access the data by xyz or ijk coordinates, 
stored as a single flat vector instead of a vector of pointers to a vector of pointers to a vector of T's. 
Contains overloaded +, -, /, *.
- Matrix and QR_fact: Contains a Matrix class with all the usual matrix operations (multiply, add and substract as 
overloaded *, +, -) plus transpose, matrix slicing by row or column, row slice operations, invert and QR factorization 
via Gram-Schmidt orthonormalization converting any invertible matrix into a product of an orthonormal Q matrix and a 
upper triangular R matrix.

