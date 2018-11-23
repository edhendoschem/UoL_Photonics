#ifndef AUX_MATHS_H_INCLUDED
#define AUX_MATHS_H_INCLUDED

#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include "matrix_fact_opt.h"

namespace Maths
{
    //Returns absolute value
    double abs_val(double const x) noexcept;
    
    //Compares equality with double values
    bool cmp_dbl(double x1, double x2, double tol = 0.00001) noexcept;
    
    //Compares the sign of two double values, returns true if they are the same
    bool sign(double const a, double const b) noexcept;
    
    //Applies the bisection method to find the root, the f_ptr must be a function that takes a
    //double as an argument
    template<typename T>
    double bisection(double const start, 
                    double const end, 
                    T f_ptr, 
                    double tol = 0.01,
                    unsigned int n_it= 10000,
                    bool pos_only = false, //Allow only positive roots
                    bool warn = false) noexcept
    {
        double x0 {start};
        double x1 {end};
        double xn {0.0};
        unsigned int count {0};
        
        while (count < n_it)
        {
            xn = {(x0 + x1) / 2.0};
            double const y0 {f_ptr(x0)};
            double const y1 {f_ptr(x1)};
            double const yn {f_ptr(xn)};
            if (warn && count % 10 == 0) std::cout<<"iteration = "<<count<<", yn = "<<yn<<'\n';
            if (Maths::abs_val(yn) < tol) 
            {
                if (pos_only) {
                    if (xn > 0.0) return xn;
                    x0 = xn;
                    x1 = end;
                    count = 0;
                } else {
                    return xn;
                }
            }
            if (Maths::sign(y0, yn))
            {
                x0 = xn;
            } else {
                x1 = xn;
            }
        
            ++count;
        
            if (count == n_it && warn)
            {
                std::cout<<"============================================\n";
                std::cout<<"Error in Maths::Bisection(), max number of iterations reached.\n";
                std::cout<<"Returning current value = "<<xn<<'\n';
                std::cout<<"Residual = "<<yn<<'\n';
                std::cout<<"============================================\n";
            }
        } //End of while loop
        
        return xn;
        
    } //End of bisection
    
    template<typename T>
    std::vector<double> bisection(std::vector<std::array<double, 2>> const& root_candidate, 
                    T f_ptr, 
                    double tol = 0.01,
                    unsigned int n_it= 10000,
                    bool pos_only = false, //Allow only positive roots
                    bool warn = false) noexcept
    {
        std::vector<double> root_values;
        
        for (auto j : root_candidate)
        {
            double const start {j[0]};
            double const end   {j[1]};
            double const x     {Maths::bisection(start, end, f_ptr, tol, n_it)};
            root_values.emplace_back(x);
        }
        
        return root_values;
    } //End of bisection vector
    
    
    //Returns the Jacobian matrix (first derivative with respect all variables for all funcs)
    template <typename T>
    Matrix<double> jacobian(std::vector<T> const& f, 
                            std::vector<double> const& x, 
                            double h = 0.00001) noexcept
    {
        std::vector<double> f_output(f.size(), 0.0);
        Matrix<double> output {f.size(), x.size()};
        
        for (auto i = 0; i < f_output.size(); ++i)
        {
            f_output[i] = f[i](x);
        }
        
        for (auto i = 0; i < f_output.size(); ++i)
        {
            for (auto j = 0; j < x.size(); ++j)
            {
                std::vector<double> x_ {x};
                x_[j] += h*x_[j];
                output(i,j) = (f[i](x_) - f_output[i]) / h;
            }
        }
        
        return output;
    } //End of jacobian function
    
    //Overload that takes matrices instead of vectors
    template <typename T>
    Matrix<double> jacobian(Matrix<T> const& f, Matrix<double> const& x, double h = 0.00001)
    {
        std::vector<T> const f_ {f.copy_vector()};
        std::vector<double> const x_ {x.copy_vector()};
        return jacobian(f_, x_, h);
    } //End of jacobian function with matrices
    
    //Overload that accepts a function that takes a vector of variables + additional 
    //non variable args (i.e. arguments that wont affect the function/derivative value)
    template <typename T, typename... ARGS>
    Matrix<double> jacobian(std::vector<T> const& f, 
                            std::vector<double> const& x, 
                            double const h = 0.0001, 
                            ARGS... args) noexcept
    {
        std::vector<double> f_output(f.size(), 0.0);
        Matrix<double> output {f.size(), x.size()};
        
        for (auto i = 0; i < f_output.size(); ++i)
        {
            f_output[i] = f[i](x, args...);
        }
        
        for (auto i = 0; i < f_output.size(); ++i)
        {
            for (auto j = 0; j < x.size(); ++j)
            {
                std::vector<double> x_ {x};
                x_[j] += h*x_[j];
                output(i,j) = (f[i](x_, args...) - f_output[i]) / h;
            }
        }
        
        return output;
    } //End of jacobian function
    
    
    template <typename T, typename... ARGS>
    Matrix<double> jacobian(Matrix<T> const& f, 
                            Matrix<double> const& x_, 
                            double const h, 
                            ARGS... args) noexcept
    {
        std::vector<double> f_output(f.size(), 0.0);
        Matrix<double> output {f.size(), x_.size()};
        std::vector<double> x {x_.copy_vector()};

        for (auto i = 0; i < f_output.size(); ++i)
        {
            f_output[i] = f[i](x, args...);
        }
        
        for (auto i = 0; i < f_output.size(); ++i)
        {
            for (auto j = 0; j < x.size(); ++j)
            {
                std::vector<double> x_ {x};
                x_[j] += h*x_[j];
                output(i,j) = (f[i](x_, args...) - f_output[i]) / h;
            }
        }
        
        return output;
    } //End of jacobian function
    
    
    
    //Jacobian Newton-Raphson, will return a vector with the first few elements as the solution
    //and the last one is the residual
    template <typename T, typename... ARGS>
    std::vector<double> jac_newton_method(std::vector<T> const& f, 
                                          std::vector<double> const& x, 
                                          double const h,
                                          double const tol,
                                          unsigned int n_it,
                                          ARGS... args) noexcept
    {
        Matrix<double> x0 {x.size(), 1, std::move(x)};
        Matrix<double> y {f.size(), 1};
        Matrix<T> f_m {f.size(), 1, std::move(f)};
        std::vector<double> output {x0.copy_vector()};
        double norm {0.0};

        
        for (auto i = 0; i < n_it; ++i)
        {
            Matrix<double> jac {Maths::jacobian(f_m.make_copy(), x0, h, args...)};
            Matrix<double> inv_jac {invert(jac)};
            
            for (auto j = 0; j < y.size(); ++j)
            {
                y(j, 0) = f_m[j](x0.copy_vector(), args...);
                norm += y(j,0) * y(j,0);
            }
            
            inv_jac * y;
            x0 - inv_jac;
            norm = sqrt(norm);
            
            if (std::isnan(norm)) 
            {
                std::cout<<"jac_newton_method Error: norm is NaN\n";
                break;
            }
            
            if (norm < tol) 
            {
                output = x0.copy_vector();
                output.push_back(norm);
                break;
            }
            
            if (i == n_it-1) 
            {
                std::cout<<"jac_newton_method warning: Max iterations reached\n";
                std::cout<<"Residual = "<<norm<<'\n';
                output = x0.copy_vector();
                output.push_back(norm);
            }
            norm = 0.0;
        }
        
        
        
        return output;
    } //End of vector jac_newton_method
    
    //Newton method with provided function pointers to the derivatives. The pointer must be
    //in the following order: {df1/dx1, df2/dx1..., dfn/dx1, df1/dx2, df2/dx2..., dfn/dx2, ...}
    template <typename T, typename... ARGS>
    std::vector<double> jac_newton_method(std::vector<T> const& f, 
                                          std::vector<double> const& x, 
                                          std::vector<T> const& df,
                                          double const tol,
                                          unsigned int n_it,
                                          ARGS... args) noexcept
    {
        Matrix<double> x0 {x.size(), 1, std::move(x)};
        Matrix<double> y {f.size(), 1};
        Matrix<T> f_m {f.size(), 1, std::move(f)};
        std::vector<double> output {x0.copy_vector()};
        double norm {0.0};

        
        for (auto i = 0; i < n_it; ++i)
        {
            Matrix<double> jac {f.size(), x.size()};
            for (auto z = 0; z < jac.size(); ++z)
            {
                jac[z] = df[z](x0.copy_vector(), args...);
            }
            
            Matrix<double> inv_jac {invert(jac)};
            
            for (auto j = 0; j < y.size(); ++j)
            {
                y(j, 0) = f_m[j](x0.copy_vector(), args...);
                norm += y(j,0) * y(j,0);
            }
            
            inv_jac * y;
            x0 - inv_jac;
            norm = sqrt(norm);
            
            if (std::isnan(norm)) 
            {
                std::cout<<"jac_newton_method Error: norm is NaN\n";
                break;
            }
            
            if (norm < tol) 
            {
                output = x0.copy_vector();
                output.push_back(norm);
                break;
            }
            
            if (i == n_it-1) 
            {
                std::cout<<"jac_newton_method warning: Max iterations reached\n";
                std::cout<<"Residual = "<<norm<<'\n';
                output = x0.copy_vector();
                output.push_back(norm);
            }
            norm = 0.0;
        }
        
        
        
        return output;
    } //End of vector jac_newton_method with provided derivatives
    
    
    //Matrix jac_newton_method
    template <typename T, typename... ARGS>
    Matrix<double> jac_newton_method(Matrix<T> const& f, 
                                          Matrix<double> const& x, 
                                          double const h,
                                          double const tol,
                                          unsigned int n_it,
                                          ARGS... args) noexcept
    {
        Matrix<double> x0 {std::move(x.make_copy())};
        Matrix<double> y {f.size(), 1};
        Matrix<T> f_m {std::move(f.make_copy())};
        Matrix<double> output {x0.size()+1, 1};
        double norm {0.0};
        
        for (auto i = 0; i < n_it; ++i)
        {
            Matrix<double> jac {Maths::jacobian(f_m.make_copy(), x0, h, args...)};
            Matrix<double> inv_jac {invert(jac)};
            for (auto j = 0; j < y.size(); ++j)
            {
                y(j, 0) = f_m[j](x0.copy_vector(), args...);
                norm += y(j,0) * y(j,0);
            }
            
            inv_jac * y;
            x0 - inv_jac;
            norm = sqrt(norm);
            
            if (std::isnan(norm)) 
            {
                std::cout<<"jac_newton_method Error: norm is NaN\n";
                break;
            }
            
            if (norm < tol) 
            {
                for (auto j = 0; j < x0.size(); ++j)
                {
                    output(j, 0) = x0(j, 0);
                }
                
                output(x0.size(), 0) = norm;
                break;
            }
            
            if (i == n_it-1) 
            {
                std::cout<<"jac_newton_method warning: Max iterations reached\n";
                std::cout<<"Residual = "<<norm<<'\n';
                for (auto j = 0; j < x0.size(); ++j)
                {
                    output(j, 0) = x0(j, 0);
                }
                
                output(x0.size(), 0) = norm;
            }
            
            norm = 0.0;
        }
        
        return output;
    } //End of jac_newton_method with matrices
    

} //namespace Maths


#endif