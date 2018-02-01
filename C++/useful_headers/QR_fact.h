#ifndef QR_FACT_H_INCLUDED
#define QR_FACT_H_INCLUDED
#include <math.h>
#include "Matrix.h"

//Auxiliary structure used to return the result from QR factorize
template<typename T>
struct QR {
    QR(Matrix<T>& Q_, Matrix<T> & R_) noexcept:
        Q{Q_}, R{R_} {};

    Matrix<T> Q; //Orthonormal matrix
    Matrix<T> R; //Upper triangular matrix
};


//Converts a matrix into the product of an orthonormal matrix (Q) and a upper triangular matrix (R). Note divisions by 0 may
//imply that the matrix is not invertible
template<typename T>
QR<T> qr_factorize(Matrix<T>& A) noexcept
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

    QR<T> output(Q, R);
    return output;
}


//Performs the specified operation between a particular matrix row and a row slice. E.G. row_slice_op(A, S, 3, '+'); adds
//to row 3 of matrix A the slice S
template<typename T>
void row_slice_op(Matrix<T>& matrix, Matrix<T>& row_slice, sz_t row, char op) noexcept
{
    if (matrix.max_cols() != row_slice.size()) {
        std::cout<<"Error in row_slice_op(): Matrix row size is different than slice size\n";
        return;
    }

    switch(op)
    {
    case '+':
        {
            for (auto j = 0; j < matrix.max_cols(); ++j) {
                matrix(row, j) = matrix(row, j) + row_slice[j];
            }
            return;
        }

    case '-':
        {
            for (auto j = 0; j < matrix.max_cols(); ++j) {
                matrix(row, j) = matrix(row, j) - row_slice[j];
            }
            return;
        }
    case '*':
        {
            for (auto j = 0; j < matrix.max_cols(); ++j) {
                matrix(row, j) = matrix(row, j) * row_slice[j];
            }
            return;
        }

    case '/':
        {
            for (auto j = 0; j < matrix.max_cols(); ++j) {
                if (!(row_slice[j]*row_slice[j] > 0 )) {
                    std::cout<<"Error in row_slice_op(): Division by zero\n";
                    return;
                }
                matrix(row, j) = matrix(row, j) / row_slice[j];
            }
            return;
        }

    case '=':
        {
            for (auto j = 0; j < matrix.max_cols(); ++j) {
                matrix(row, j) = row_slice[j];
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


//inverts the target matrix. Note: Division by 0 errors imply that the matrix may not be invertible
template<typename T>
Matrix<T> invert(Matrix<T> const& A_) noexcept
{
    //Check if it is a square matrix
    if (A_.max_cols()!=A_.max_rows()) {
        std::cout<<"Error in invert(): Not a square matrix (m != n)\n";
        Matrix<T> output;
        return output;
    }

    Matrix<T> output = A_; //Create a copy of the matrix to convert to identity
    make_identity(output); //Convert it to identity matrix
    Matrix<T> A = A_; //Create another copy of the original matrix

    //Make A a lower triangular matrix and its diagonal 1
    for (auto i = 0; i < A.max_rows(); ++i) {
        //create a pivot
        T const mult = A(i,i);
        Matrix<T> slice_A = A.slice(i, 'r');
        Matrix<T> slice_output = output.slice(i, 'r');
        slice_A = slice_A / mult;
        slice_output = slice_output / mult;
        row_slice_op(A, slice_A, i, '=');
        row_slice_op(output, slice_output, i, '=');

        //Make all numbers below pivot zero
        for (auto j = i+1; j < A.max_cols(); ++j) {
            if (i == A.max_rows()-1) break;
            Matrix<T> temp_slice_1 = slice_A * A(j, i);
            Matrix<T> temp_slice_2 = slice_output * A(j, i);
            row_slice_op(A, temp_slice_1, j, '-');
            row_slice_op(output, temp_slice_2, j, '-');
        }
    }

    sz_t m_ = A.max_rows()-1;
    sz_t n_ = A.max_cols()-1;
    std::cout<<A<<'\n';
    //Make the upper triangular part of the matrix 0 and obtain the inverse in output. Start from the bottom
    for (auto i = 0; i < A.max_rows(); ++i) {

        Matrix<T> slice_A = A.slice(m_-i, 'r');
        Matrix<T> slice_output = output.slice(m_-i, 'r');

        //Use bottom pivot and make elements above it zero
        for (auto j = i+1; j < A.max_cols(); ++j) {
            if (i == A.max_rows()-1) break; //Ignore first row (i.e. last value of the loop)
            Matrix<T> temp_slice_1 = slice_A * A(n_-j, m_-i);
            Matrix<T> temp_slice_2 = slice_output * A(n_-j, m_-i);
            row_slice_op(A, temp_slice_1, n_-j, '-');
            row_slice_op(output, temp_slice_2, n_-j, '-');
        }
    }

    return output;
}


#endif // QR_FACT_H_INCLUDED
