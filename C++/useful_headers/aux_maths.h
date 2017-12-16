#ifndef AUX_MATHS_H_INCLUDED
#define AUX_MATHS_H_INCLUDED

#include <cmath> //Allows the use of cmath functions
#include <complex> //Allows the use of complex numbers and functions
#include <numeric>
#include <vector> //enables std::vector
#include <array> //enables std::array
#include <random> //Enables random numbers
#include <iostream>
#include "Flat_vec.h" //Enables Flat_vec class


//Mersenne Twister recommended initialization
inline std::mersenne_twister_engine<std::uint_fast64_t, 64, 312, 156, 31,
                             0xb5026f5aa96619e9, 29,
                             0x5555555555555555, 17,
                             0x71d67fffeda60000, 37,
                             0xfff7eee000000000, 43, 6364136223846793005> Mersenne_Twister{};
using Rand_e = decltype(Mersenne_Twister);


//Absolute value that returns double
inline double abs_d(double val) {
    if (val < 0.0) {
        return -1.0 * val;
    } else {
        return val;
    }
}

//Integrates functions using Monte Carlo, res = number of points, n_cycles is number of cycles
//Free function calculator
template<typename T>
decltype(auto) integrate_fn(std::vector<double> a, std::vector<double> b,  T& fn_ptr,
                            unsigned long long res = 800000, unsigned long long n_cycles = 1,
                            Rand_e rand_eng = Mersenne_Twister)
                            noexcept {
    std::vector<std::uniform_real_distribution<>> distributions;
    std::vector<double> args(a.size(), 0.0);
    double multiplier = 1.0;
    double res_d = static_cast<double>(res);
    double n_cycles_d = static_cast<double>(n_cycles);
    double sum = 0;
    double sum_total = 0;

    for (auto i = 0; i < a.size(); ++i) {
        std::uniform_real_distribution<> unif_real{a[i], b[i]};
        distributions.push_back(unif_real);
        multiplier *= (b[i] - a[i]);
    }


    for (auto cycles = 0; cycles < n_cycles; ++cycles) {
        for (unsigned long long i = 0; i < res; ++i) {
            for (auto j = 0; j < distributions.size(); ++j) {
                double const arg = distributions[j](rand_eng);
                args[j] = arg;
            }

            sum += fn_ptr(args);
        }

        sum /= res_d;
        sum_total += sum;
        sum = 0;
    }

    sum_total = sum_total * multiplier / n_cycles_d;
    return sum_total;
}


//Static member function calculator
template<typename T, typename R>
decltype(auto) integrate_fn(std::vector<double> a, std::vector<double> b,  T& fn_ptr, R& obj_ptr,
                            unsigned long long res = 800000, unsigned long long n_cycles = 1,
                            Rand_e rand_eng = Mersenne_Twister)
                            noexcept {
    std::vector<std::uniform_real_distribution<>> distributions;
    std::vector<double> args(a.size(), 0.0);
    double multiplier = 1.0;
    double res_d = static_cast<double>(res);
    double n_cycles_d = static_cast<double>(n_cycles);
    double sum = 0;
    double sum_total = 0;

    for (auto i = 0; i < a.size(); ++i) {
        std::uniform_real_distribution<> unif_real{a[i], b[i]};
        distributions.push_back(unif_real);
        multiplier *= (b[i] - a[i]);
    }


    for (auto cycles = 0; cycles < n_cycles; ++cycles) {
        for (unsigned long long i = 0; i < res; ++i) {
            for (auto j = 0; j < distributions.size(); ++j) {
                double const arg = distributions[j](rand_eng);
                args[j] = arg;
            }
            //std::cout<<"fn_ptr val = "<<fn_ptr(args, obj_ptr)<<'\n';
            sum += fn_ptr(args, obj_ptr);
        }

        sum /= res_d;
        sum_total += sum;
        sum = 0;
    }

    sum_total = sum_total * multiplier / n_cycles_d;
    return sum_total;
}


//Returns the result of ordinary differential equations dy/dx = f(x,y) with y(x_0) = y0

//Runge-Kutta for free functions
template<typename T>
double runge_kutta(std::vector<double> params, T& fn_ptr, unsigned long long h = 100000) noexcept {
    double const h_d = static_cast<double> (h);
    double const x0 = params[0];
    double const y0 = params[1];
    double const x_max = params[2];
    double h_sz = x_max / h_d;
    double x{x0}, y{y0}, dxy {0}, k1{0}, k2{0}, k3{0}, k4{0};

    for (auto i = 0; i < h; ++i) {
        dxy += fn_ptr(x, y);
        k1 = h_sz * fn_ptr(x, y);
        k2 = h_sz * fn_ptr(x + 0.5 * h_sz, y + 0.5*k1);
        k3 = h_sz * fn_ptr(x + 0.5 * h_sz, y + 0.5*k2);
        k4 = h_sz * fn_ptr(x + h_sz, y + k3);
        y = y + (k1/6.0) + (k2/3.0) + (k3/3.0) + (k4/6.0);
        x = x + h_sz;
    }

    return y;
}


//Runge-Kutta for static functions inside classes/structs
template<typename T, typename R>
double runge_kutta(std::vector<double> params, T& fn_ptr, R& obj_ptr, unsigned long long h = 100000) noexcept {
    double const h_d = static_cast<double> (h);
    double const x0 = params[0];
    double const y0 = params[1];
    double const x_max = params[2];
    double h_sz = x_max / h_d;
    double x{x0}, y{y0}, dxy {0}, k1{0}, k2{0}, k3{0}, k4{0};

    for (auto i = 0; i < h; ++i) {
        dxy += fn_ptr(x, y, obj_ptr);
        k1 = h_sz * fn_ptr(x, y, obj_ptr);
        k2 = h_sz * fn_ptr(x + 0.5 * h_sz, y + 0.5*k1, obj_ptr);
        k3 = h_sz * fn_ptr(x + 0.5 * h_sz, y + 0.5*k2, obj_ptr);
        k4 = h_sz * fn_ptr(x + h_sz, y + k3, obj_ptr);
        y = y + (k1/6.0) + (k2/3.0) + (k3/3.0) + (k4/6.0);
        x = x + h_sz;
    }

    return y;
}
#endif // AUX_MATHS_H_INCLUDED

//Bisection method

struct Res_type { //Auxiliary struct to signal whether a result is in the valid state
    Res_type(double _value, bool _valid)
    : value{_value}, valid{_valid} {};

    double value {0.0};
    bool valid {false};
};

template<typename T, typename... ARGS>
Res_type bisection(std::vector<double> x_lims, T& fn_ptr,
                 double tol, unsigned long long steps, ARGS&... args) noexcept {

    double a_a = x_lims[0];
    double b_b = x_lims[1];
    double c_c = 0.0;
    double f_a_a = fn_ptr(a_a, args...);
    double f_b_b = fn_ptr(b_b, args...);
    double f_c_c = 0.0;
    if (abs_d(f_a_a) < tol) {
        Res_type result {a_a, true};
        return result;
    }
    if (abs_d(f_b_b) < tol) {
        Res_type result {b_b, true};
        return result;
    }

    if (f_a_a > 0 && f_b_b > 0 || f_a_a < 0 && f_b_b < 0) {
        std::cout<<"Bisection(): Upper and lower intervals have the same sign, no root finding possible, returning 0.0\n";
        Res_type result {0.0, false};
        return result;
    }

    for (auto i = 0; i < steps; ++i) {
        c_c = (a_a + b_b) / 2.0;
        f_c_c = fn_ptr(c_c, (args)...);

        if (abs_d(f_c_c) < tol) {
            Res_type result {c_c, true};
            return result;
        }
        if (f_c_c < 0) {
            a_a = c_c;
        } else {
            b_b = c_c;
        }
    }

    Res_type result {c_c, false};
    return result;
}
