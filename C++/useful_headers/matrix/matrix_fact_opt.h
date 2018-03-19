#ifndef MATRIX_FACT_OPT_H_INCLUDED
#define MATRIX_FACT_OPT_H_INCLUDED
#include <math.h>
#include "matrix_opt.h"

//Cholesky factorization
//Auxiliary function to determine symmetry
template<typename T>
bool is_symmetric(Matrix<T> const& m_sys_, sz_t limit = 1000000) noexcept
{
    auto is_sym = [] (sz_t starts, sz_t ends, Matrix<T> const& m_sys)
    {

        if (m_sys.max_cols()  == m_sys.max_rows())
        {
            for (auto i = 0; i < m_sys.max_rows(); ++i) {
                for (auto j = i+1; j < m_sys.max_cols(); ++j) {
                    if (m_sys(i,j)  == m_sys(j,i)) {continue;}
                    else {return false;}
                }
            }

        } else {
            return false;
        }

        return true;
    };

    if (limit < m_sys_.size() && limit > THREADS) {
        std::vector<std::future<bool>> future_res;
        sz_t new_start = 0;
        sz_t new_end = 0;
        for (auto i = 0; i < THREADS; ++i) {
            new_end =(i+1) * m_sys_.max_rows()/THREADS;
            std::cout<<"new_start = "<<new_start<<"new_end = "<<new_end<<'\n';
            std::future<bool> res = std::async(std::launch::async, is_sym, new_start, new_end, std::ref(m_sys_));
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
        return is_sym(0, m_sys_.max_rows(), m_sys_);
    }
}


//Cholesky factorization, output matrix pair containing a lower triangular matrix and an upper triangular
template<typename T>
Matrix_pair<T> cholesky_factorize_opt(Matrix<T> const& A_, sz_t limit = 1000000) noexcept
{
    static std::atomic<std::size_t> curr_num {4};
    if (is_symmetric(A_) == false) {
        Matrix_pair<T> result {};
        std::cout<<"Error in cholesky_factorize(): Matrix not symmetric. Returning default Matrix_pair\n";
        return result;
    }

    Matrix<T> B {A_.max_rows(), A_.max_cols()}; //Lower triangular

    auto chol_fac = [] (sz_t starts, sz_t ends, sz_t i, Matrix<T> const& A, Matrix<T>& B)
    {
        std::cout<<"starts = "<<starts<<", ends = "<<ends<<'\n';
        T sum_result = 0.0;
        ends = ends + i > A.max_cols() ? A.max_cols() : ends + i;

        for (auto j = starts+i; j < ends; ++j) {
            for (auto k = 0; k < i; ++k) {
                sum_result += B(i,k) * B(i,k);
            }

            B(i,i) = sqrt(A(i,i) - sum_result);
            sum_result = 0;

            for (auto k = 0; k < i; ++k) {
                sum_result += B(i,k) * B(j,k);
            }

            B(j,i) = (1.0/B(i,i)) * (A(i,j) - sum_result);
            sum_result = 0;
        }

        return;
    };

    if (limit < A_.size() && limit > THREADS) {
        for (auto i = 0; i < A_.max_rows(); ++i) {
            launch_task(chol_fac, 0, A_.max_cols(), 0, i, std::ref(A_), std::ref(B));
        }
    } else {
        for (auto i = 0; i < A_.max_rows(); ++i) {
            chol_fac(0, A_.max_cols(), i, A_, B);
        }
    }

    Matrix<T> C {A_.max_rows(), A_.max_cols()}; //Upper triangular
    C = std::move(B.make_copy());
    transpose(C);

    Matrix_pair<T> result {std::move(B), std::move(C)};

    return result;
}


//Modified version of Matrix.h's row_slice_op to be used only with lu_factorize function
template<typename T>
void row_slice_op2(Matrix_pair<T>& m_sys, Matrix_pair<T>& row_slice, sz_t row, sz_t curr_row, char op) noexcept
{
    if (m_sys.first.max_cols()  != row_slice.first.size()  ||
        m_sys.second.max_cols() != row_slice.second.size())
    {
        std::cout<<"Error in row_slice_op(): Matrix(ces) row size(s) is/are different than slice size(s)\n";
        return;
    }

    switch(op)
    {
    case '+':
        {

            for (auto j = 0; j < m_sys.second.max_cols(); ++j) {
                m_sys.second(row, j) = m_sys.second(row, j) + row_slice.second[j];
            }

            for (auto k = 0; k < curr_row; ++k) {
                row_slice.first[k] = 0.0;
            }

            for (auto j = 0; j < m_sys.first.max_cols(); ++j) {
                m_sys.first(row, j) = m_sys.first(row, j) - row_slice.first[j];
            }

            return;
        }

    case '-':
        {

            for (auto j = 0; j < m_sys.second.max_cols(); ++j) {
                m_sys.second(row, j) = m_sys.second(row, j) - row_slice.second[j];

            }

            for (auto k = 0; k < curr_row; ++k) {
                row_slice.first[k] = 0.0;
            }

            for (auto j = 0; j < m_sys.first.max_cols(); ++j) {

                m_sys.first(row, j) = m_sys.first(row, j) + row_slice.first[j];

            }

            return;
        }

    default:
        {
            std::cout<<"Error in row_slice_op(): Unknown operation\n";
            return;
        }
    }

    return;
}


//LU factorize, outputs a matrix pair containing a lower triangular and a upper triangular matrix
template<typename T>
Matrix_pair<T> lu_factorize(Matrix<T> const& A, sz_t limit = 1000000) noexcept
{

    Matrix<T> L(A.max_rows(), A.max_cols());
    make_identity(L);
    Matrix<T> U {std::move(A.make_copy())};
    Matrix_pair<T> output {std::move(L), std::move(U)};
    sz_t m = output.first.max_rows();
    sz_t n = output.first.max_cols();

    auto lu_fac = [] (sz_t starts, sz_t ends, sz_t j, Matrix_pair<T>& output)
    {
        ends = ends + j > output.first.max_rows() ? output.first.max_rows() : ends + j;

        for (auto i = starts+j; i < ends; ++i) {
            Matrix_pair<T> row_slices = output.slice(j, 'r');
            T const mult = output.second(i,j) / row_slices.second[j];
            mult * row_slices;
            row_slice_op2(output, row_slices, i, j, '-');
        }

        return;
    };

    if (limit < output.first.size() && limit > THREADS) {
        for (auto j = 0; j < n; ++j) {
            launch_task(lu_fac, 1, output.first.max_rows(), 0, j, std::ref(output));
        }
    } else {
        for (auto j = 0; j < n; ++j) {
            lu_fac(1, output.first.max_rows(), j, output);
        }
    }

    return output;
}

//Calculate the determinant value by using LU factorization
template<typename T>
T determinant(Matrix_pair<T> const& LU) noexcept
{
    T output = 1;

    for (auto i = 0; i < LU.first.max_cols(); ++i) {
        output *= LU.second(i,i);
    }

    return output;
}


//Calculate the determinant value by using LU factorization
template<typename T>
T determinant(Matrix_system<T> const& LU) noexcept
{
    T output = 1;

    for (auto i = 0; i < LU.first.max_cols(); ++i) {
        output *= LU.second(i,i);
    }

    return output;
}



//Converts a matrix into the product of an orthonormal matrix (Q) and a upper triangular matrix (R). Note divisions by 0 may
//imply that the matrix is not invertible

//Aux function for qr_factorize
template<typename T>
T find_col_modulus(Matrix<T> const& M, size_t const k) {
    T result {0};

    for (auto i = 0; i < M.max_rows(); ++i) {
        result += (M(i, k) * M(i, k));
    }

    result = sqrt(result);
    return result;
}


template<typename T>
Matrix_pair<T> qr_factorize(Matrix<T> const& A, sz_t limit = 1000000) noexcept
{
    //Auxiliary vector of column matrices required for the Gram-Schmidt orthonormalization process
    Matrix<T> u {A.max_rows(), A.max_cols()};
    Matrix<T> e {A.max_rows(), A.max_cols()};

    //Calculate the values for e and u matrices
    T temp = find_col_modulus(A, 0);
    Matrix<T> aj {A.slice(0, 'c')};
    col_slice_op(u, aj, 0, '=');
    aj / temp;
    col_slice_op(e, aj, 0, '=');

    for (auto j = 0; j < A.max_cols(); ++j) {
        Matrix<T> aj {A.slice(j, 'c')};
        Matrix<T> aj_orig {A.slice(j, 'c')};

        for (auto k = 0; k < j; ++k) {
            Matrix<T> ek {e.slice(k, 'c')};
            Matrix<T> ek_t {e.slice(k, 'c')};
            Matrix<T> aj_temp {aj_orig.make_copy()};
            transpose(aj_temp);
            aj_temp * ek_t;
            ek * aj_temp;
            aj - ek;
        }

        col_slice_op(u, aj, j, '=');
        T temp = find_col_modulus(u, j);
        Matrix<T> ej {u.slice(j, 'c')};
        ej / temp;
        col_slice_op(e, ej, j, '=');
    }

    Matrix<T> Q {std::move(e)};
    Matrix<T> R {A.max_rows(), A.max_cols()};

    for (auto i = 0; i < A.max_rows(); ++i) {
        Matrix<T> e_temp {std::move(Q.slice(i, 'c'))};

        for (auto j = i; j < A.max_cols(); ++j) {
            Matrix<T> A_temp {std::move(A.slice(j, 'c'))};
            transpose(A_temp, limit);
            A_temp * e_temp;
            R(i,j) = A_temp(0,0);
        }
    }

    Matrix_pair<T> output {std::move(Q), std::move(R)};
    return output;
}


//Solvers
//Use LU factorization with. The output is a Matrix containing a single column with the solutions
template<typename T>
Matrix<T> lu_solve(Matrix<T> const& M, Matrix<T> const& C, sz_t limit = 1000000) noexcept
{
    if (M.max_rows() != C.max_rows()) {
        std::cout<<"Error in lu_solve: constants matrix' number of rows differs from coefficient matrix' number of rows.";
        std::cout<<" Returning default constructed matrix\n";
        Matrix<T> output {1,1};
        return output;
    }

    Matrix<T> solution {std::move(C.make_copy())};
    Matrix_pair<T> output {std::move(lu_factorize(M, limit))};

    for (auto i = 0; i < C.size(); ++i) {
        for (auto j = 0; j < i; ++j) {
            solution[i] -= solution[j] * output.first(i,j);
        }
    }


    for (auto i = static_cast<long long> (C.size()-1); i >= 0; --i) {
        T temp = solution[i];
        for (auto j = static_cast<long long> (C.size()-1); j > i; --j) {
            temp -= solution[j] * output.second(i,j);
        }

        temp = temp / output.second(i,i);
        solution[i] = temp;
    }

    return solution;
}


//Use Cholesky to solve an equation system. The output is a Matrix containing a single column with the solutions
template<typename T>
Matrix<T> cholesky_solve(Matrix<T> const& M, Matrix<T> const& C, sz_t limit = 1000000) noexcept
{
    if (M.max_rows() != C.max_rows()) {
        std::cout<<"Error in cholesky_solve: constants matrix' number of rows differs from coefficient matrix' number";
        std::cout<<"of rows. Returning default constructed matrix\n";
        Matrix<T> output {1,1};
        return output;
    }

    Matrix<T> solution {std::move(C.make_copy())};
    Matrix_pair<T> output {std::move(lu_factorize(M))};

    for (auto i = 0; i < C.size(); ++i) {
        T temp = solution[i];

        for (auto j = 0; j < i; ++j) {
            temp -= solution[j] * output.first(i,j);
        }

        temp = temp / output.first(i, i);
        solution[i] = temp;
    }


    for (auto i = static_cast<long long> (C.size()-1); i >= 0; --i) {
        T temp = solution[i];

        for (auto j = static_cast<long long> (C.size()-1); j > i; --j) {
            temp -= solution[j] * output.second(i,j);
        }

        temp = temp / output.second(i,i);
        solution[i] = temp;
    }

    return solution;
}


template<typename T>
Matrix<T> qr_solve(Matrix<T> const& M, Matrix<T> const& C, sz_t limit = 1000000) noexcept
{
    if (M.max_rows() != C.max_rows()) {
        std::cout<<"Error in qr_solve: constants matrix' number of rows differs from coefficient matrix' number of rows.";
        std::cout<<" Returning default constructed matrix\n";
        Matrix<T> output {1,1};
        return output;
    }

    Matrix<T> solution {std::move(C.make_copy())};
    Matrix_pair<T> output {std::move(lu_factorize(M))};

    for (auto i = 0; i < C.size(); ++i) {
        T temp = solution[i];

        for (auto j = 0; j < i; ++j) {
            temp -= solution[j] * output.first(i,j);
        }

        temp = temp / output.first(i, i);
        solution[i] = temp;
    }

    for (auto i = static_cast<long long> (C.size()-1); i >= 0; --i) {
        T temp = solution[i];

        for (auto j = static_cast<long long> (C.size()-1); j > i; --j) {
            temp -= solution[j] * output.second(i,j);
        }

        temp = temp / output.second(i,i);
        solution[i] = temp;
    }

    return solution;
}
#endif // MATRIX_FACT_OPT_H_INCLUDED
