#ifndef MATRIX_FACT_H_INCLUDED
#define MATRIX_FACT_H_INCLUDED

#include "matrix.h"
#include <algorithm>

namespace Maths
{

    //Exchanges row with target row within the same matrix
    template <typename T>
    void swap_rows(Matrix<T>& A, sz_t const row, sz_t const target_row, sz_t const limit = 500000) noexcept
    {
        auto max_c {A.max_cols()};
        auto inner_launcher = [&] (sz_t const st, sz_t ed)
        {
            for (auto j = st; j < ed; ++j)
            {
                std::swap(A(row, j), A(target_row, j));
            }
        };
        
        if (A.size() > limit && max_c > THREADS)
        {
            launch_task(inner_launcher, static_cast<sz_t>(0), max_c, 0);
        } else {
            inner_launcher(static_cast<sz_t>(0), max_c);
        }
        
        return;
    }

    //Exchanges row with target row within the same matrix for all matrices
    template <typename T>
    void swap_rows(Matrix_pair<T>& A, sz_t const row, sz_t const target_row, sz_t const limit = 500000) noexcept
    {
        swap_rows(A.first, row, target_row, limit);
        swap_rows(A.second, row, target_row, limit);
        
        return;
    }


    //Exchanges row with target row within the same matrix for all matrices
    template <typename T>
    void swap_rows(Matrix_system<T>& A, sz_t const row, sz_t const target_row, sz_t const limit = 500000) noexcept
    {
        swap_rows(A.first, row, target_row, limit);
        swap_rows(A.second, row, target_row, limit);
        swap_rows(A.third, row, target_row, limit);
        
        return;
    }

    
    //Helper function to automatically find the largest pivot in a particular column 
    //and return its row index. Note it begins searching at row = col as it is assumed the
    //previous max rows would have been switched
    template <typename T>
    sz_t find_pivot_row(Matrix<T> const&A, sz_t const col, sz_t const limit = 500000) noexcept
    {
        auto max_r {A.max_rows()};
        
        auto inner_launcher = [&] (sz_t const st, sz_t const ed)
        {
            T largest;
            sz_t idx{col};
            
            for (auto i = st ; i < ed; ++i)
            {
                if (i == st)
                {
                    largest = std::abs(A(i, col));
                    idx = i;
                    continue;
                    
                }

                if(std::abs(A(i, col)) > largest)
                {
                    largest = std::abs(A(i, col));
                    idx = i;
                }
            }
            
            return idx;
        };
        
        if (A.size() > limit && max_r > THREADS)
        {
            std::vector<std::future<sz_t>> future_res;
            std::vector<std::array<sz_t, 2>> intervals {split_interval(col, max_r)};
            
            for (auto i = 0; i < THREADS; ++i)
            {
                std::future<sz_t> res = std::async(std::launch::async, 
                                                    inner_launcher, 
                                                    intervals[i][0], 
                                                    intervals[i][1]);
                future_res.push_back(std::move(res));    
            }
            
            sz_t output {col};
            
            for (auto i = 0; i < THREADS; ++i)
            {
                sz_t const temp {future_res[i].get()};
                if (std::abs(A(temp, col)) > std::abs(A(output, col))) output = temp;
            }
            
            return output;
        } else {
            sz_t output {inner_launcher(col, max_r)};
            return output;
        }
    } //End of find_pivot_row
    

    //Cholesky factorization
    //Auxiliary function to determine symmetry
    template<typename T>
    bool is_symmetric(Matrix<T> const& m_sys_, sz_t const limit = 500000) noexcept
    {
        if (m_sys_.max_cols() != m_sys_.max_rows()) return false;
    
        auto inner_launcher = [&m_sys_] (sz_t starts, sz_t ends)
        {
            auto max_r = m_sys_.max_rows();
            auto max_c = m_sys_.max_cols();
        
            if (max_r  == max_c)
            {
                for (auto i = 0; i < max_r; ++i) {
                    for (auto j = i+1; j < max_c; ++j) {
                        if (m_sys_(i,j)  == m_sys_(j,i)) {continue;}
                        else {return false;}
                    }
                }
            } else {
                return false;
            }

            return true;
        };

        if (m_sys_.size() > limit && limit > THREADS) {
            std::vector<std::future<bool>> future_res;
            
            sz_t new_start = 0;
            sz_t new_end = 0;
            std::vector<std::array<sz_t, 2>> intervals {split_interval(static_cast<sz_t>(0), m_sys_.max_rows())};
        
            for (auto i = 0; i < THREADS; ++i) {
                std::future<bool> res = std::async(std::launch::async, 
                                                    inner_launcher, 
                                                    intervals[i][0], 
                                                    intervals[i][1]);
                future_res.push_back(std::move(res));
                new_start = new_end;
            }

            std::vector<bool> results_(future_res.size());
            for (auto i = 0; i < future_res.size(); ++i) {
                results_[i] = future_res[i].get();
            }

            for (bool elem : results_) {
                if (elem ==  false) return false;
            }

            return true;

        } else {
            return inner_launcher(0, m_sys_.max_rows());
        }
    } //end of is_symmetric


    //Cholesky factorization, output matrix pair containing a lower triangular matrix and an upper triangular
    template<typename T>
    Matrix_pair<T> cholesky_factorize(Matrix<T> const& A_, sz_t const limit = 500000)
    {
        if (is_symmetric(A_) == false) {
            std::string const e{"Error in cholesky_factorize(): Matrix not symmetric.\n"};
            throw Error{e, Error_type::matrix_not_symmetric};
        }
        
        try
        {
            auto max_r {A_.max_rows()};
            auto max_c {A_.max_cols()};
    
            Matrix<T> B {max_r, max_c}; //Lower triangular

            auto inner_launcher = [&] (sz_t const st, sz_t ed, sz_t const i)
            {
                T sum_result = 0.0;
                ed = ed + i > max_c ? max_c : ed + i;

                for (auto j = st+i; j < ed; ++j) {
                    for (auto k = 0; k < i; ++k) {
                        sum_result += B(i,k) * B(i,k);
                    }

                    B(i,i) = sqrt(A_(i,i) - sum_result);
                    sum_result = 0;

                    for (auto k = 0; k < i; ++k) {
                        sum_result += B(i,k) * B(j,k);
                    }

                    B(j,i) = (1.0/B(i,i)) * (A_(i,j) - sum_result);
                    sum_result = 0;
                }

                return;
            };

            if (A_.size() >limit && limit > THREADS) {
                for (auto i = 0; i < max_r; ++i) {
                    launch_task(inner_launcher, 0, max_c, 0, i);
                }
            } else {
                for (auto i = 0; i < max_r; ++i) {
                    inner_launcher(0, max_c, i);
                }
            }

            Matrix<T> C {max_r, max_c}; //Upper triangular
            C = B;
            C = transpose(C);

            Matrix_pair<T> result {std::move(B), std::move(C)};

            return result;
        }
        catch (Error& e)
        {
            std::string const prev_e {e.what()};
            std::string const this_e {"Error in cholesky_factorize() -> "+prev_e};
            throw Error{this_e, e.et};
        }
    } //End of Cholesky


    //LU factorize, outputs a matrix pair containing a lower triangular and a upper triangular 
    //matrix. PA = LU where P is a permutation matrix and A the original matrix
    template<typename T>
    Matrix_pair<T> lu_factorize(Matrix<T> const& A, sz_t limit = 500000)
    {
        try
        {
            auto max_r {A.max_rows()};
            auto max_c {A.max_cols()};
        
            Matrix<T> L {max_r, max_c};
            Matrix<T> U {A};
            Matrix<T> permutation {make_identity<T>(max_r, max_c)};
            Matrix_pair<T> output {std::move(L), std::move(U)};
            T const neg {static_cast<T>(-1.0)};
        
            for (auto i = 0; i < max_r-1; ++i)
            {
                sz_t const pivot_idx {find_pivot_row(output.second, i, limit)};
                if (pivot_idx != i) 
                {
                    swap_rows(output, i, pivot_idx, limit);
                    swap_rows(permutation, i, pivot_idx, limit);
                }
            
                T const curr_pivot {output.second(i,i)};
            
                for (auto j = i+1; j < max_c; ++j)
                {
                    Matrix<T> curr_slice {output.second.slice(i, 'r')};
                    T const mult {neg * output.second(j,i) / curr_pivot};
                    curr_slice = curr_slice * mult;
                    output.second = row_slice_op(output.second, curr_slice, j, '+');
                    output.first(j,i) = neg * mult;
                }
            }
        

            for (auto i = 0; i < max_r; ++i)
            {
                output.first(i,i) = 1.0;
            }

            return output;
        }
        catch (Error& e)
        {
            std::string const prev_e {e.what()};
            std::string const this_e {"Error in lu_factorize() -> "+prev_e};
            throw Error{this_e, e.et};
        }
    } //End of lu factorize


    //Helper options enum
    enum class Opt {Cholesky, LU};
    //Calculate the determinant value by using LU factorization
    template<typename T>
    T determinant(Matrix<T> const& matrix, Opt const o = Opt::LU, sz_t const limit = 500000)
    {
        try 
        {
            switch (o)
            {
                case (Opt::LU):
                {
                    T output = static_cast<T>(1);
                    Matrix_pair<T> lu {lu_factorize(matrix)};
                    auto max_r {lu.second.max_rows()};
    
                    for (auto i = 0; i < max_r; ++i) {
                        output *= lu.second(i,i);
                    }

                    return output;
                }
                
                case (Opt::Cholesky):
                {
                    T output = static_cast<T>(1);
                    Matrix_pair<T> chol {cholesky_factorize(matrix)};
                    auto max_r {chol.second.max_rows()};
    
                    for (auto i = 0; i < max_r; ++i) {
                        output *= chol.second(i,i);
                    }
                    
                    output *= output;

                    return output;
                }
                
                default:
                {
                    std::string const e {"Error in determinant(): Unknown option\n"};
                    throw Error{e, Error_type::invalid_option};
                }
            }
            
        }
        catch(Error& e)
        {
            std::string const prev_e {e.what()};
            std::string const this_e {"Error in determinant() -> "+prev_e};
            throw Error{this_e, e.et};
        }
    } //end of determinant




    //Converts a matrix into the product of an orthonormal matrix (Q) and a upper triangular matrix (R). Note divisions by 0 may
    //imply that the matrix is not invertible

    //Aux function for qr_factorize
    template<typename T>
    T find_col_modulus(Matrix<T> const& M, size_t const k) 
    {
        try
        {
            T result {static_cast<T>(0)};

            for (auto i = 0; i < M.max_rows(); ++i) {
                result += (M(i, k) * M(i, k));
            }

            result = sqrt(result);
            return result;
        }
        catch(Error& e)
        {
            std::string const prev_e {e.what()};
            std::string const curr_e {"Error in find_col_modulus() -> "+prev_e};
            throw Error{curr_e, e.et};
        }
    }


    template<typename T>
    Matrix_pair<T> qr_factorize(Matrix<T> const& A, sz_t limit = 500000)
    {
        try
        {
            auto max_r {A.max_rows()};
            auto max_c {A.max_cols()};
            //Auxiliary vector of column matrices required for the Gram-Schmidt 
            //orthonormalization procedure
            Matrix<T> u {max_r, max_c};
            Matrix<T> e {max_r, max_c};
            
            //Calculate the values for e and u matrices
            T temp = find_col_modulus(A, 0);
 
            Matrix<T> aj {A.slice(0, 'c')};
            u  = col_slice_op(u, aj, 0, '=');
            aj = aj / temp;
            e  = col_slice_op(e, aj, 0, '=');
            
            for (auto j = 0; j < max_c; ++j) {
                Matrix<T> aj {A.slice(j, 'c')};
                Matrix<T> aj_orig {A.slice(j, 'c')};

                for (auto k = 0; k < j; ++k) {
                    Matrix<T> ek {e.slice(k, 'c')};
                    Matrix<T> ek_t {e.slice(k, 'c')};
                    Matrix<T> aj_temp {aj_orig};
                    aj_temp = transpose(aj_temp);
                    aj_temp = aj_temp * ek_t;
                    ek = ek * aj_temp;
                    aj = aj - ek;
                }

                u = col_slice_op(u, aj, j, '=');
                T temp = find_col_modulus(u, j);
                Matrix<T> ej {u.slice(j, 'c')};
                ej = ej / temp;
                e = col_slice_op(e, ej, j, '=');
            }

            Matrix<T> Q {std::move(e)};
            Matrix<T> R {max_r, max_c};

            for (auto i = 0; i < max_r; ++i) {
                Matrix<T> e_temp {Q.slice(i, 'c')};

                for (auto j = i; j < max_c; ++j) {
                    Matrix<T> A_temp {A.slice(j, 'c')};
                    A_temp = transpose(A_temp, limit);
                    A_temp = A_temp * e_temp;
                    R(i,j) = A_temp(0,0);
                }
            }

            Matrix_pair<T> output {std::move(Q), std::move(R)};
            return output;
        }
        catch(Error& e)
        {
            std::string const prev_e {e.what()};
            std::string const curr_e {"Error in qr_factorize() -> "+prev_e};
            throw Error{curr_e, e.et};
        }
    } //End of qr_factorize


    //Finds matrix eigenvalues and outputs a sorted vector with the results
    template<typename T>
    std::vector<T> find_eigenvalues(Matrix<T> const& matrix, sz_t const iterations = 10) noexcept
    {
        try
        {
            if (matrix.max_cols() != matrix.max_rows())
            {
                std::string const e {"Error in find_eigenvalues(): Non square matrix\n"};
                throw Error{e, Error_type::matrix_not_square};
            }
            
            std::vector<T> output;
            Matrix<T> m {matrix};
            Matrix_pair<T> mp; 
            auto max_r {m.max_rows()};
            
            for (auto i = 0;  i < iterations; ++i)
            {
                mp = qr_factorize(m);
                m = mp.second * mp.first;
            }
            
            for (auto i = 0; i < max_r; ++i)
            {
                output.emplace_back(m(i,i));
            }
            
            std::sort(output.begin(), output.end());
            return output;
        }
        catch(Error& e)
        {
            std::string const prev_e {e.what()};
            std::string const curr_e {"Error in find_eigenvalues() -> " + prev_e};
            throw Error{curr_e, e.et};
        }
    } //End of find_eigenvalues


    //inverts the target matrix. Note: Division by 0 errors imply that the matrix may not be invertible
    template<typename T>
    Matrix<T> invert(Matrix<T> const& matrix, sz_t const limit = 500000)
    {
        T det {determinant<T>(matrix)};
        if (is_zero<T>(det))
        {
            std::string const e {"Error in invert(): Zero value determinant\n"};
            throw Error{e, Error_type::zero_value_determinant};
        }
    
        try
        {
            auto max_r {matrix.max_rows()};
            auto max_c {matrix.max_cols()};
            
            Matrix_pair<T> mp {matrix, make_identity<T>(max_r, max_c)};
            T const neg {static_cast<T>(-1.0)};
            
            for (auto i = 0; i < max_r; ++i)
            {
                sz_t const pivot_idx {find_pivot_row(mp.first, i, limit)};
                if (pivot_idx != i) 
                {
                    swap_rows(mp, i, pivot_idx, limit);
                }
                
                T const curr_pivot {mp.first(i,i)};
                
                for (auto j = 0; j < max_c; ++j)
                {
                    mp.first(i,j)  = mp.first(i,j) / curr_pivot;
                    mp.second(i,j) = mp.second(i,j) / curr_pivot;
                }
                
                //Break once it reaches to the bottom
                if (i == max_r - 1) break;
                
                //Eliminating all lower triangular values
                for (auto j = i+1; j < max_c; ++j)
                {
                    Matrix_pair<T> curr_slice {mp.slice(i, 'r')};
                    T const mult {neg * mp.first(j,i)};
                    curr_slice * mult;
                    row_slice_op(mp, curr_slice, j, '+');
                }
            }
            
            //Eliminating upper triangular values
            for (auto i = 1; i < max_r; ++i)
            {
                for (auto j = 0; j < i; ++j)
                {
                    T const mult {mp.first(j, i)};
                    Matrix_pair<T> curr_slice {mp.slice(i, 'r')};
                    curr_slice * (neg * mult);
                    row_slice_op(mp, curr_slice, j, '+');
                }
            }
            
            Matrix<T> output {std::move(mp.second)};
            return output;
        }
        catch(Error& e)
        {
            std::string const prev_e {e.what()};
            std::string const curr_e {"Error in invert() -> " + prev_e};
            throw Error{curr_e, e.et};
        }
    } //End of matrix inversion

} //End of namespace Maths
#endif