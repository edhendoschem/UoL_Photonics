#include "aux_maths.h"

//Returns absolute value
double Maths::abs_val(double const x) noexcept
{
    if (x < 0.0) return -x;
    return x;
}
    
//Compares equality with double values
bool Maths::cmp_dbl(double x1, double x2, double tol) noexcept
{
    if ((Maths::abs_val((x1-x2))) < tol) return true;
    return false;
}
    
//Compares the sign of two double values, returns true if they are the same
bool Maths::sign(double const a, double const b) noexcept
{
    if ((a >= 0.0 && b >= 0.0) || (a < 0.0 && b < 0.0)) return true;
    return false;
}
