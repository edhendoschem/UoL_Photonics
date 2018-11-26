#ifndef AUX_MATHS_H_INCLUDED
#define AUX_MATHS_H_INCLUDED
#include "matrix_fact.h"
#include <random>

namespace Maths
{
    class Random_engine
    {
    public:
        
        Random_engine()
        {
            std::random_device rd;
            std::array<double, std::mt19937::state_size> seed_data{};
            std::generate(std::begin(seed_data), std::end(seed_data), std::ref(rd));
            std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
            engine.seed(seq);
        }
        
        std::mt19937& operator () () noexcept
        {
            return engine;
        }
        
    private:
        std::mt19937 engine;
    }; //Random_engine
    

    template <typename T, typename F_ptr, typename... ARGS>
    std::vector<T> derivative(F_ptr const& f, std::vector<T> const& x, double const h, ARGS... args)
    {
        if (is_zero(h))
        {
            std::string const e {"Error in derivative(): differential h cannot be zero\n"};
            throw Error{e};
        }
        
        sz_t const x_size {x.size()};
        std::vector<T> output(x_size, 0.0);
        T const f_x {f(x, args...)};
        T const h_t {static_cast<T>(h)};
        
        for (auto i = 0; i < x_size; ++i)
        {
            std::vector<T> x_mod {x};
            T x_h;
            if (is_zero(h_t * x_mod[i]))
            {
                x_h = h_t;
            } 
            else 
            {
                x_h = h_t * x_mod[i];
            }
            
            x_mod[i] += x_h;
            T const f_x_h {f(x_mod, args...)};
            output[i] = (f_x_h - f_x) / x_h;
        }
        
        return output;
        
    } //End of derivaive vector version


    template <typename T, typename F_ptr, typename... ARGS>
    Matrix<T> derivative(F_ptr const& f, Matrix<T> const& x, double const h, ARGS... args)
    {
        std::vector<T> result {derivative(f, x.get_data(), h, args...)};
        Matrix<T> output {static_cast<sz_t>(1), x.size(), std::move(result)};
        return output;
    }
    
    //Matrix variable
    template <typename T, typename F_ptr, typename... ARGS>
    Matrix<T> jacobian(std::vector<F_ptr> const& f, Matrix<T> const& x, double const h, ARGS... args)
    {
        auto max_r {f.size()};
        auto max_c {x.size()};
        
        try
        {
            Matrix<T> output{max_r, max_c};
            sz_t count {0};
        
            for (auto const& i : f)
            {
                Matrix<T> const row {std::move(derivative(i, x, h, args...))};
                output = row_slice_op(output, row, count, '=');
                ++count;
            }
        
            return output;
        }
        catch(Error& e)
        {
            std::string const prev_e {e.what()};
            std::string const curr_e {"Error in jacobian -> " + prev_e};
            throw Error{e};
        }
    }
    
    //Vector variables
    template <typename T, typename F_ptr, typename... ARGS>
    Matrix<T> jacobian(std::vector<F_ptr> const& f, std::vector<T> const& x, double const h, ARGS... args)
    {
        try
        {
            Matrix<T> x_mat {x.size(), 1, x};
            Matrix<T> output {jacobian(f, x_mat, h, args...)};
            
            return output;
        }
        catch(Error& e)
        {
            std::string const prev_e {e.what()};
            std::string const curr_e {"Error in jacobian -> " + prev_e};
            throw Error{e};
        }
    }    
    
    
    
    
}//End of namespace Maths




#endif