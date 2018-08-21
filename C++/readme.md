# Table of contents

## 1. [new_simulation](#it1)
## 2. [useful_headers](#it2)
## 3. [EDWA_sim] (#it3)

## new_simulation <a name = "it1"></a>
This is a finite difference time domain method solver for Maxwell's Equations, featuring second order absorbing boundary
in 1D and 2D grids, also the ability to easily add Perfect electrical conductors (PEC) and dielectric structures and 
be able to choose from a variety of sources.

<b>Note:</b> the text files beginning with plot_script_  are GNUPLOT scripts used to quickly visualize the data
produced by the simulation

## useful_headers <a name = "it2"></a>
Contains:

- aux_maths: Contains Bisection method, Jacobian matrix calculation, Newton method, other useful functions
- Flat_vec: A vector with defined access operator in order to be able to access the data by xyz or ijk coordinates, 
stored as a single flat vector instead of a vector of pointers to a vector of pointers to a vector of T's. 
Contains overloaded +, -, /, *.
- matrix_opt and matrix_fact_opt Contains: 
	* A move only Matrix class with all the usual matrix operations (multiply, add and substract as 
overloaded *, +, - operators) which are stored in the first matrix for memory efficiency e.g a+b will add b to a, the result
will be stored in a. A "make_copy()" member function is provided should the need for copying arise
	* Other operations such as transpose, matrix slicing by row or column, row slice operations, column slice operations and invert
	* QR factorization via Gram-Schmidt orthonormalization converting any invertible matrix into a product of an orthonormal Q matrix 
and a upper triangular R matrix 
	* LU factorization which converts a given matrix into the product of a lower triangular with an upper triangular and matrix determinant 
	* Cholesky factorize which is a highly robust and more numerically efficient LU decomposition, only applicable for symmetric, positive definite matrices. 
<b>Note</b>: If the size of the matrix is 
larger than 1000x1000, parallel processing will be used in the following functions: lu_factorize, cholesky_factorize, transpose
	* Matrix linear equation system solvers using the aforementioned methods lu_solve, cholesky_solve and qr_solve

## 2018-06-13_EDWA_sim (in progress) <a name = "it3"></a>
Four level Erbium, 2 level Ytterbium simulation which takes into account up-conversion, cross-relaxation between Ytterbium and Erbium. Currently in progress. Requires
matrix and aux_maths