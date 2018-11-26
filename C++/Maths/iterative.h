#ifndef ITERATIVE_H_INCLUDED
#define ITERATIVE_H_INCLUDED
#include "aux_maths.h"

namespace Maths
{
    template <typename T>
    struct Iter_params
    {
        Iter_params() : h {1e-6}, tol {1e-7}, n_it {1000} {}
        Iter_params(T const h_, T const tol_, sz_t const n_it_):
            h {h_}, tol {tol_}, n_it {n_it_} {}
        T h;
        T tol;
        sz_t n_it;
    };
    
    
    //Applies newton method to find the root of an equation system made of functions that take n 
    //variables and returns a single variable each Vector version
    template <typename T, typename F_ptr, typename... ARGS>
    std::vector<T> newton_method (std::vector<F_ptr> const& f, std::vector<T> const& x, Iter_params<T> const& p, ARGS... args)
    {
        try 
        {
            T const h {p.h};
            T const tol {p.tol};
            sz_t const n_it {p.n_it};
            
            std::vector<T> output {x};
            Matrix<T> x_curr {x.size(), 1, x};
            
            //Special case for overconstrained problems
            if (f.size() > x.size())
            {
                for (auto n = 0; n < n_it; ++n)
                {
                    Matrix<T> f_curr {f.size(), 1};
                    sz_t count {0};
                    T total {0};
                
                    for (auto k : f)
                    {
                        f_curr[count] = f[count](x_curr.get_data(), args...);
                        total += std::abs(f_curr[count]);
                        ++count;
                    
                    }
                
                    if (total < tol) break;
                    
                    Matrix<T> const jac{jacobian(f, x_curr, h, args...)};
                    Matrix<T> const t_jac {transpose(jac)};
                    Matrix<T> const intermediate {t_jac * jac};
                    Matrix<T> const inv_jac {invert(intermediate)};
                    Matrix<T> const pseudo_inv {inv_jac * t_jac};
                    
                    x_curr = x_curr - (pseudo_inv * f_curr);
                }

                output = x_curr.get_data();
                return output;
            
            } //End special case for overconstrained problems
            
            //Regular case (Equal number of functions)
            for (auto n = 0; n < n_it; ++n)
            {
                Matrix<T> f_curr {f.size(), 1};
                sz_t count {0};
                T total {0};
                
                for (auto k : f)
                {
                    f_curr[count] = f[count](x_curr.get_data(), args...);
                    total += std::abs(f_curr[count]);
                    ++count;
                    
                }
                
                if (total < tol) break;
                    
                Matrix<T> const jac{jacobian(f, x_curr, h, args...)};
                Matrix<T> const inv_jac {invert(jac)};
                x_curr = x_curr - (inv_jac * f_curr);
            }

            output = x_curr.get_data();
            return output;
        }
        catch (Error& e)
        {
            std::string const prev_e {e.what()};
            std::string const curr_e {"Error in newton_method() -> " + prev_e};
            throw Error{curr_e};
        }
        
    } //End of newton_method
    
    
    //Applies newton method to find the root of an equation system made of functions that take n 
    //variables and returns a single variable each Matrix version
    template <typename T, typename F_ptr, typename... ARGS>
    Matrix<T> newton_method (std::vector<F_ptr> const& f, Matrix<T> const& x, Iter_params<T> const& p, ARGS... args)
    {
        std::vector<T> const x_vec {x}
        std::vector<T> out_vec {newton_method(f, x_vec, p, args...)};
        Matrix<T> const output {x.size(), 1, std::move(out_vec)};
        return output;
    }
} //End of namespace maths

#endif