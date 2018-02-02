#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED
#include <iostream>
#include <vector>

using sz_t = std::size_t;
using sz_t = std::size_t;


template<typename T>
class Matrix
{
private:
    std::vector<T> data; //Container of matrix data
    sz_t m; //Number of rows
    sz_t n; //Number of columns

public:
    //Constructors
    //Default constructor
    Matrix() noexcept:
        m{1}, n{1}
    {
        data.resize(1, 0.0);
    }


    //Construct matrix given dimensions
    Matrix(sz_t m_, sz_t n_) noexcept:
        m{m_}, n{n_}
    {
        if (m*n > 0) {
            data.resize(m*n, 0.0);
        } else {
            set_cols(1);
            set_rows(1);
        }
    }


    //Copy constructor
    Matrix(Matrix const& matrix) noexcept:
        m{matrix.max_rows()}, n{matrix.max_cols()}
    {
        data = std::move(matrix.get_data());
    }


    //Overloaded operators
    //Access an element by inputting the row and column value, constat version
    T const& operator () (sz_t const i, sz_t const j) const noexcept
    {
        if (i < m && j < n) {
            return data[i + (m * j)];
        } else {
            if (i < m) {
                std::cout<<"Error in operator(): row index "<<i<<" exceeds m-1. Returning first element\n";
            } else {
                std::cout<<"Error in operator(): column index "<<j<<" exceeds n-1. Returning first element\n";
            }
            return data[0];
        }
    }


    //Access an element by inputting the row and column value, non-const version
    T& operator () (sz_t const i, sz_t const j) noexcept
    {
        if (i < m && j < n) {
            return data[i + (m * j)];
        } else {
            if (i < m) {
                std::cout<<"Error in operator(): row index "<<i<<" exceeds m-1. Returning first element\n";
            } else {
                std::cout<<"Error in operator(): column index "<<j<<" exceeds n-1. Returning first element\n";
            }
            return data[0];
        }
    }


    //Access the vector that contains the elements directly by index const version
    T const& operator [] (sz_t const i) const noexcept
    {
        if (i < data.size()) {
            return data[i];
        } else {
            std::cout<<"Error in operator[]: Index "<<i<<" exceeds maximum size. Returning first element\n";
            return data[0];
        }
    }


    //Access the vector that contains the elements directly by index non-const version
    T& operator [] (sz_t const i) noexcept
    {
        if (i < data.size()) {
            return data[i];
        } else {
            std::cout<<"Error in operator[]: Index "<<i<<" exceeds maximum size. Returning first element\n";
            return data[0];
        }
    }


    //Move-copy operator
    Matrix<T>& operator = (Matrix<T> const& matrix) = default;


    //Utility functions
    //Return a column slice or row slice. Note: This slice is a copy and not part of the matrix
    Matrix<T> slice (sz_t const k, char c) const noexcept
    {
        switch(c)
        {
        case 'c':
            {
                Matrix<T> output(m,1);
                for (auto i = 0; i < m; ++i) {
                    output[i] = data[i + k*m];
                }

                return output;
            }

        case 'r':
            {
                Matrix<T> output(1,n);
                for (auto j = 0; j < n; ++j) {
                    output[j] = data[k + j*m];
                }

                return output;
            }

        default:
            {
                std::cout<<"Error in operator(): invalid option. returning default constructed matrix\n";
                Matrix<T> output{};
                return output;
            }

        }
    }


    //Return the size of the underlying vector container
    std::size_t size() const noexcept
    {
        return data.size();
    }


    //Change the underliying vector
    void set_data(std::vector<T> const& data_) noexcept
    {
        data = data_;
    }


    //Return the underlying vector
    std::vector<T> const& get_data() const noexcept
    {
        return data;
    }


    //Returns the number of rows in the matrix
    sz_t max_rows() const noexcept
    {
        return m;
    }


    //Returns the number of columns in the matrix
    sz_t max_cols() const noexcept
    {
        return n;
    }


    //Changes the number of rows in the vector, either adding zeros at the end (increase) or truncating the matrix (decrease)
    void set_rows(sz_t m_) noexcept
    {
        if (m == m_) return;
        if (m_ == 0) {
            std::cout<<"Error in set_rows(): Number of rows cannot be 0\n";
            return;
        }

        sz_t difference = m > m_ ? m-m_:m_-m;
        std::vector<T> new_vector;
        new_vector.resize(m_*n, 0);
        if (m_ > m) {
            auto k = 0;
            auto counter = m-1;
            for (auto i = 0; i < data.size(); ++i) {
                new_vector[k] = data[i];
                if(i == counter) {
                    k += 1+difference;
                    counter += m;
                } else {
                    ++k;
                }
            }
        } else {
            auto k = 0;
            auto counter = m-1;
            for (auto i = 0; i < data.size(); ++i) {
                new_vector[k] = data[i];
                if(i == counter) {
                    k += 1-difference;
                    counter += m;
                } else {
                    ++k;
                }
            }
        }

        m = m_;
        data = std::move(new_vector);
    }


    //Changes the number of columns in the vector, either adding zeros at the end (increase) or truncating the matrix (decrease)
    void set_cols(sz_t n_) noexcept
    {
        if (n == n_) return;
        if (n_ == 0) {
            std::cout<<"Error in set_cols(): Number of columns cannot be 0\n";
            return;
        }

        sz_t difference = n > n_ ? n-n_:n_-n;
        std::vector<T> new_vector;
        new_vector.resize(m*n_, 0);

        if (n_ > n) {
            for (auto j = 0; j < data.size(); ++j) {
                new_vector[j] = data[j];
            }
        } else {
            for (auto j = 0; j < data.size() - difference; ++j) {
                new_vector[j] = data[j];
            }
        }

        n = n_;
        data = std::move(new_vector);
    }

};


//Other operators
//Overloading << operator to help print the matrix
template<typename T>
std::ostream& operator << (std::ostream& os, Matrix<T> matrix) noexcept
{
    for (auto i = 0; i < matrix.max_rows(); ++i) {
        for (auto j = 0; j < matrix.max_cols(); ++j) {
            if (j != matrix.max_cols()-1) {
                os<<matrix(i,j)<<'\t'<<'\t';
            } else {
                os<<matrix(i,j)<<'\n';
            }
        }
    }

    return os;
}


//Specialized << operator, to print a vector of column matrixes
template<typename T>
std::ostream& operator << (std::ostream& os, std::vector<Matrix<T>> vec) noexcept
{

    for (auto j = 0; j < vec[0].size(); ++j) {
        for (auto i = 0; i < vec.size(); ++i) {
            if (i != vec.size() - 1) {
                os<<vec[i][j]<<'\t';
            } else {
                os<<vec[i][j]<<'\n';
            }
        }
    }

    return os;
}


//Adds a constant to every element of the matrix, outputs a new matrix
template<typename T>
Matrix<T> operator + (Matrix<T> const& matrix_a, T const c) noexcept
{
    Matrix<T> output = matrix_a;
    for (auto i = 0; i < output.max_rows(); ++i) {
        for (auto j = 0; j < output.max_cols(); ++j) {
            output(i,j) += c;
        }
    }

    return output;
}


//Adds a constant to every element of the matrix, outputs a new matrix
template<typename T>
Matrix<T> operator + (T const c, Matrix<T> const& matrix_a) noexcept
{
    Matrix<T> output = matrix_a;
    for (auto i = 0; i < output.max_rows(); ++i) {
        for (auto j = 0; j < output.max_cols(); ++j) {
            output(i,j) += c;
        }
    }

    return output;
}


//Adds two matrices together
template<typename T>
Matrix<T> operator + (Matrix<T> const& matrix_a, Matrix<T> const& matrix_b) noexcept
{
    if (matrix_a.max_rows() != matrix_b.max_rows() || matrix_a.max_cols() != matrix_b.max_cols()) {
        std::cout<<"Error in operator +: Incompatible matrix dimensions. Returning default constructed matrix\n";
        Matrix<T> output{};
        return output;
    }

    Matrix<T> output(matrix_a.max_rows(), matrix_a.max_cols());

    for (auto i = 0; i < output.max_rows(); ++i) {
        for (auto j = 0; j < output.max_cols(); ++j) {
            output(i,j) = matrix_a(i,j) + matrix_b(i,j);
        }
    }

    return output;
}


//Substracts a constant to all elements of the matrix, the output is a new matrix
template<typename T>
Matrix<T> operator - (Matrix<T> const& matrix_a, T const c) noexcept
{
    Matrix<T> output = matrix_a;
    for (auto i = 0; i < output.max_rows(); ++i) {
        for (auto j = 0; j < output.max_cols(); ++j) {
            output(i,j) -= c;
        }
    }

    return output;
}


//Substracts a constant to all elements of the matrix, the output is a new matrix
template<typename T>
Matrix<T> operator - (T const c, Matrix<T> const& matrix_a) noexcept
{
    Matrix<T> output = matrix_a;
    for (auto i = 0; i < output.max_rows(); ++i) {
        for (auto j = 0; j < output.max_cols(); ++j) {
            output(i,j) = c - output(i,j);
        }
    }

    return output;
}


//Substracts two matrices
template<typename T>
Matrix<T> operator - (Matrix<T> const& matrix_a, Matrix<T> const& matrix_b) noexcept
{
    if (matrix_a.max_rows() != matrix_b.max_rows() || matrix_a.max_cols() != matrix_b.max_cols()) {
        std::cout<<"Error in operator +: Incompatible matrix dimensions. Returning default constructed matrix\n";
        Matrix<T> output;
        return output;
    }

    Matrix<T> output(matrix_a.max_rows(), matrix_a.max_cols());

    for (auto i = 0; i < output.max_rows(); ++i) {
        for (auto j = 0; j < output.max_cols(); ++j) {
            output(i,j) = matrix_a(i,j) - matrix_b(i,j);
        }
    }

    return output;
}


//Multiplies two matrices of compatible dimensions
template<typename T>
Matrix<T> operator * (Matrix<T> const& matrix_a, Matrix<T> const& matrix_b) noexcept
{
    if (matrix_a.max_cols() != matrix_b.max_rows()) {
        std::cout<<"Error in operator *: Incompatible matrix dimensions. Returning default constructed matrix\n";
        Matrix<T> output;
        return output;
    }

    sz_t m_ = matrix_a.max_rows();
    sz_t n_ = matrix_b.max_cols();
    Matrix<T> output(m_, n_);

    for (auto i = 0; i < matrix_a.max_rows(); ++i) {
        for (auto j = 0; j < matrix_b.max_cols(); ++j) {
            for (auto k = 0; k < matrix_a.max_cols(); ++k) {
                output(i, j) += matrix_a(i, k) * matrix_b(k, j);
            }
        }
    }

    return output;
}


//Multiplies all elements of the matrix by a constant
template<typename T>
Matrix<T> operator * (Matrix<T> const& matrix_a, T const c) noexcept
{
    Matrix<T> output{matrix_a};

    for (auto i = 0; i < output.max_rows(); ++i) {
        for (auto j = 0; j < output.max_cols(); ++j) {
            output(i,j) *= c;
        }
    }

    return output;
}


//Multiplies all elements of the matrix by a constant
template<typename T>
Matrix<T> operator * (T const c, Matrix<T> const& matrix_a) noexcept
{
    Matrix<T> output{matrix_a};

    for (auto i = 0; i < output.max_rows(); ++i) {
        for (auto j = 0; j < output.max_cols(); ++j) {
            output(i,j) *= c;
        }
    }

    return output;
}


//Divides all elements of the matrix by a constant
template<typename T>
Matrix<T> operator / (Matrix<T> const& matrix_a, T const c) noexcept
{
    if (!(c*c > 0)) {
        std::cout<<"Error in operator /: Division by zero encountered. Returning default constructed matrix\n";
        Matrix<T> output{};
        return output;
    }

    Matrix<T> output{matrix_a};

    for (auto i = 0; i < output.max_rows(); ++i) {
        for (auto j = 0; j < output.max_cols(); ++j) {
            output(i,j) /= c;
        }
    }

    return output;
}


//Utility functions
//Switches the rows with the columns of the matrix
template<typename T>
Matrix<T> transpose(Matrix<T> const& matrix_a) noexcept
{
    Matrix<T> output(matrix_a.max_cols(), matrix_a.max_rows());

    for (auto i = 0; i < matrix_a.max_rows(); ++i) {
        for (auto j = 0; j < matrix_a.max_cols(); ++j) {
            output(j,i) = matrix_a(i,j);
        }
    }

    return output;
}


//Creates a matrix with 0 in all the elements except those in the diagonal wich are equal to 1
template<typename T>
void make_identity(Matrix<T>& matrix_a) noexcept
{
    for (auto i = 0; i < matrix_a.max_rows(); ++i) {
        for (auto j = 0; j < matrix_a.max_cols(); ++j) {
            if (i == j) {
                matrix_a(i,j) = 1;
            } else {
                matrix_a(i,j) = 0;
            }
        }
    }
}

#endif // MATRIX_H_INCLUDED
