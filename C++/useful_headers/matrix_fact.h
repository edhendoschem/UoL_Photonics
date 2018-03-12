#ifndef MATRIX_FACT_H_INCLUDED
#define MATRIX_FACT_H_INCLUDED
#include <math.h>
#include "Matrix.h"

//Cholesky factorization
//Auxiliary function to determine symmetry
template<typename T>
bool is_symmetric(Matrix<T> const& m_sys) noexcept
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
}


//Cholesky factorization, output matrix pair
template<typename T>
Matrix_pair<T> cholesky_factorize(Matrix<T> const& A) noexcept
{
    if (is_symmetric(A) == false) {
        Matrix_pair<T> result {};
        std::cout<<"Error in cholesky_factorize(): Matrix not symmetric. Returning default Matrix_pair\n";
        return result;
    }

    Matrix<T> B {A.max_rows(), A.max_cols()}; //Lower triangular
    Matrix<T> C {A.max_rows(), A.max_cols()}; //Upper triangular
    T sum_result = 0.0;

    for (auto i = 0; i < A.max_rows(); ++i) {
        for (auto j = i; j < A.max_cols(); ++j) {
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
    }

    C = transpose(B);
    Matrix_pair<T> result {B, C};

    return result;
}


//Modified version of Matrix.h's row_slice_op to be used only with lu_factorize function
template<typename T>
void row_slice_op2(Matrix_system<T>& m_sys, Matrix_system<T>& row_slice, sz_t row, sz_t curr_row, char op) noexcept
{
    if (m_sys.first.max_cols()  != row_slice.first.size()  ||
        m_sys.second.max_cols() != row_slice.second.size() ||
        m_sys.third.max_cols()  != row_slice.third.size())
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

            for (auto j = 0; j < m_sys.third.max_cols(); ++j) {
                m_sys.third(row, j) = m_sys.third(row, j) + row_slice.third[j];

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

            for (auto j = 0; j < m_sys.third.max_cols(); ++j) {
                m_sys.third(row, j) = m_sys.third(row, j) - row_slice.third[j];

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


//Lower Upper (LU) factorization. Converts a matrix A into the product of a lower triangular matrix and an upper triangular matrix
template<typename T>
Matrix_pair<T> lu_factorize(Matrix<T>& A) noexcept
{

    Matrix<T> L(A.max_rows(), A.max_cols());
    make_identity(L);
    Matrix<T> U = A;
    Matrix_pair<T> output {L, U};

    sz_t m = output.first.max_rows();
    sz_t n = output.first.max_cols();

    for (auto j = 0; j < n; ++j) {
        for (auto i = j+1; i < m; ++i) {
            Matrix_pair<T> row_slices = output.slice(j, 'r');
            T const mult = output.second(i,j) / row_slices.second[j];
            row_slices = mult * row_slices;
            row_slice_op2(output, row_slices, i, j, '-');
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


//Converts a matrix into the product of an orthonormal matrix (Q) and a upper triangular matrix (R). Note divisions by 0 may
//imply that the matrix is not invertible
template<typename T>
Matrix_pair<T> qr_factorize(Matrix<T>& A) noexcept
{
    Matrix<T> Q(A.max_rows(), A.max_cols());
    Matrix<T> R = Q;

    //Auxiliary vector of column matrices required for the Gram-Schmidt orthonormalization process
    std::vector<Matrix<T>> a;
    std::vector<Matrix<T>> e;
    std::vector<Matrix<T>> u;

    //Fill a vector with column matrices
    for (auto j = 0; j< A.max_cols(); ++j) {
        Matrix<T> col_vec_a(A.max_rows(), 1);
        for(auto i = 0; i < A.max_rows(); ++i) {
            col_vec_a(i,0) = A(i,j);
        }
        a.push_back(col_vec_a);
    }

    e = a; //starting values same as a
    u = a; //starting values same as a

    //Calculate the values for e and u matrices
    for (auto j = 0; j< A.max_cols(); ++j) {
        T mod = 0;

        if (j == 0) {
            for(auto i = 0; i < A.max_rows(); ++i) {
                mod += (u[j][i] * u[j][i]);
            }

            mod = sqrt(mod);
            e[j] = u[j]/mod;
        } else {
            mod = 0;

            u[j] = a[j];
            for (auto k = 0; k < j; ++k) {
                Matrix<T> temporary = transpose(a[j]) * e[k];
                T mult = temporary[0];
                u[j] = u[j] - (mult * e[k]);
            }

            for(auto i = 0; i < A.max_rows(); ++i) {
                mod += (u[j][i] * u[j][i]);
            }

            mod = sqrt(mod);
            e[j] = u[j] / mod;

        }
    }

    //Complete the orthonormal matrix
    for (auto i = 0; i < e.size(); ++i) {
        for (auto j = 0; j < e[0].size(); ++j) {
            Q(j,i) = e[i][j];
        }
    }

    //Complete the upper triangular matrix
    for (auto i = 0; i < R.max_rows(); ++i) {
        for (auto j = 0; j < R.max_cols(); ++j) {
            if (j >= i) {
                Matrix<T> temp = transpose(A.slice(j,'c')) * Q.slice(i,'c');
                T member = temp[0];
                R(i,j) = member;
            } else {
                R(i,j) = 0;
            }
        }
    }

    Matrix_pair<T> output(Q, R);
    return output;
}



#endif // MATRIX_FACT_H_INCLUDED
