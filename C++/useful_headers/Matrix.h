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


    //Copy assignment operator
    Matrix<T>& operator = (Matrix<T> const& matrix) = default;


    //Utility functions
    //Return a column slice or row slice. Note: This slice is a copy and not part of the matrix
    Matrix<T> slice (sz_t const k, char c) const noexcept
    {
        switch(c)
        {
        case 'c':
            {


                if ((m-1)+k*m < data.size()) {
                    Matrix<T> output(m,1);
                    for (auto i = 0; i < m; ++i) {
                        output[i] = data[i+k*m];
                    }

                    return output;

                } else {
                    std::cout<<"Error in .slice(): Invalid column. Returning default constructed matrix\n";
                    Matrix<T> output{};
                    return output;
                }
            }

        case 'r':
            {
                if (k+(n-1)*m < data.size()) {
                    Matrix<T> output(1,n);
                    for (auto j = 0; j < n; ++j) {
                        output[j] = data[k + j*m];
                    }

                    return output;
                } else {
                    std::cout<<"Error in .slice(): Invalid row. Returning default row matrix\n";
                    Matrix<T> output{};
                    return output;
                }
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


    //Change the underlying vector
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
            auto counter = m_-1;
            for (auto i = 0; i < new_vector.size(); ++i) {
                new_vector[i] = data[k];
                if(i == counter) {
                    k += 1+difference;
                    counter += m_;
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

        std::vector<T> new_vector;
        new_vector.resize(m*n_, 0);

        if (n_ > n) {
            for (auto j = 0; j < data.size(); ++j) {
                new_vector[j] = data[j];
            }
        } else {
            for (auto j = 0; j < new_vector.size(); ++j) {
                new_vector[j] = data[j];
            }
        }

        n = n_;
        data = std::move(new_vector);
    }


    //Vector iterators
    auto cbegin() const noexcept
    {
        return data.cbegin();
    }


    auto begin() noexcept
    {
        return data.begin();
    }


    auto cend() const noexcept
    {
        return data.cend();
    }


    auto end() noexcept
    {
        return data.end();
    }

};


//Helper structs
//Auxiliary struct for returning pair of matrices
template<typename T>
struct Matrix_pair {
    Matrix<T> first;
    Matrix<T> second;

    //Constructors
    //Default constructor
    Matrix_pair() noexcept:
        first{Matrix<T>()}, second{Matrix<T>()} {};


    //Constructor from two matrices
    Matrix_pair(Matrix<T> const& first_, Matrix<T> const& second_) noexcept:
        first{first_}, second{second_} {};


    //Copy constructor
    Matrix_pair(Matrix_pair<T> const& S_) noexcept:
        first{S_.first}, second{S_.second} {};


    //Copy assignment operator
    Matrix_pair<T>& operator = (Matrix_pair<T> const& B) = default;


    //Returns a column slice or row slice of the entire system. Note: This slice is a copy and not part of the matrix
    Matrix_pair<T> slice (sz_t const k, char c) const noexcept
    {

        switch(c)
        {
        case 'c':
            {
                Matrix<T> A = first.slice(k, 'c');
                Matrix<T> B = second.slice(k, 'c');
                Matrix_pair<T> output {A, B};
                return output;
            }

        case 'r':
            {
                Matrix<T> A = first.slice(k, 'r');
                Matrix<T> B = second.slice(k, 'r');
                Matrix_pair<T> output {A, B};
                return output;
            }

        default:
            {
                std::cout<<"Error in operator(): invalid option. returning default constructed matrix pair\n";
                Matrix_pair<T> output{};
                return output;
            }
        }
    }
};


//Auxiliary struct to represent system of matrices
template<typename T>
struct Matrix_system {
    Matrix<T> first;
    Matrix<T> second;
    Matrix<T> third;


    //Constructors
    //Default constructor
    Matrix_system() noexcept:
        first{Matrix<T>()}, second{Matrix<T>()}, third{Matrix<T>()} {};

    //Constructor for homogeneous matrices system
    Matrix_system(Matrix<T> const& first_, Matrix<T> const& second_) noexcept:
        first{first_}, second{second_}, third{Matrix<T>(second_.max_rows(), 1)} {};

    //Constructor for three matrices system
    Matrix_system(Matrix<T> const& first_, Matrix<T> const& second_, Matrix<T> const& third_) noexcept:
        first{first_}, second{second_}, third{third_} {};

    //Copy constructor
    Matrix_system(Matrix_system<T> const& S_) noexcept:
        first{S_.first}, second{S_.second}, third{S_.third} {};


    //Copy assignment operator
    Matrix_system<T>& operator = (Matrix_system<T> const& B) = default;


    //Returns a column slice or row slice of the entire system. Note: This slice is a copy and not part of the matrix
    Matrix_system<T> slice (sz_t const k, char c) const noexcept
    {

        switch(c)
        {
        case 'c':
            {
                Matrix<T> A = first.slice(k, 'c');
                Matrix<T> B = second.slice(k, 'c');
                Matrix<T> C = third.slice(k, 'c');
                Matrix_system<T> output {A, B, C};
                return output;
            }

        case 'r':
            {
                Matrix<T> A = first.slice(k, 'r');
                Matrix<T> B = second.slice(k, 'r');
                Matrix<T> C = third.slice(k, 'r');
                Matrix_system<T> output {A, B, C};
                return output;
            }

        default:
            {
                std::cout<<"Error in operator(): invalid option. returning default constructed matrix system\n";
                Matrix_system<T> output{};
                return output;
            }

        }
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
                os<<matrix(i,j)<<' '<<' ';
            } else {
                os<<matrix(i,j)<<'\n';
            }
        }
    }

    return os;
}


//Overloaded << operator to print a matrix pair
template<typename T>
std::ostream& operator << (std::ostream& os, Matrix_pair<T> matrix) noexcept
{
    std::cout<<"First matrix\n";
    for (auto i = 0; i < matrix.first.max_rows(); ++i) {
        for (auto j = 0; j < matrix.first.max_cols(); ++j) {
            if (j != matrix.first.max_cols()-1) {
                os<<matrix.first(i,j)<<' '<<' ';
            } else {
                os<<matrix.first(i,j)<<'\n';
            }
        }
    }

    std::cout<<"Second matrix\n";
    for (auto i = 0; i < matrix.second.max_rows(); ++i) {
        for (auto j = 0; j < matrix.second.max_cols(); ++j) {
            if (j != matrix.second.max_cols()-1) {
                os<<matrix.second(i,j)<<' '<<' ';
            } else {
                os<<matrix.second(i,j)<<'\n';
            }
        }
    }

    return os;
}


//Overloaded << operator to print a matrix system
template<typename T>
std::ostream& operator << (std::ostream& os, Matrix_system<T> matrix) noexcept
{
    std::cout<<"First matrix\n";
    for (auto i = 0; i < matrix.first.max_rows(); ++i) {
        for (auto j = 0; j < matrix.first.max_cols(); ++j) {
            if (j != matrix.first.max_cols()-1) {
                os<<matrix.first(i,j)<<' '<<' ';
            } else {
                os<<matrix.first(i,j)<<'\n';
            }
        }
    }

    std::cout<<"Second matrix\n";
    for (auto i = 0; i < matrix.second.max_rows(); ++i) {
        for (auto j = 0; j < matrix.second.max_cols(); ++j) {
            if (j != matrix.second.max_cols()-1) {
                os<<matrix.second(i,j)<<' '<<' ';
            } else {
                os<<matrix.second(i,j)<<'\n';
            }
        }
    }

    std::cout<<"Third matrix\n";
    for (auto i = 0; i < matrix.third.max_rows(); ++i) {
        for (auto j = 0; j < matrix.third.max_cols(); ++j) {
            if (j != matrix.third.max_cols()-1) {
                os<<matrix.third(i,j)<<' '<<' ';
            } else {
                os<<matrix.third(i,j)<<'\n';
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
    if (matrix_b.max_cols() == 1 && matrix_b.max_rows() == 1) {
        Matrix<T> output {matrix_a};
        T const temp = matrix_b(0,0);
        output = output + temp;
        return output;
    }

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
    if (matrix_b.max_cols() == 1 && matrix_b.max_rows() == 1) {
        Matrix<T> output {matrix_a};
        T const temp = matrix_b(0,0);
        output = output - temp;
        return output;
    }

    if (matrix_a.max_rows() != matrix_b.max_rows() || matrix_a.max_cols() != matrix_b.max_cols()) {
        std::cout<<"Error in operator -: Incompatible matrix dimensions. Returning default constructed matrix\n";
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


//Multiplies two matrices of compatible dimensions
template<typename T>
Matrix<T> operator * (Matrix<T> const& matrix_a, Matrix<T> const& matrix_b) noexcept
{
    if (matrix_b.max_cols() == 1 && matrix_b.max_rows() == 1) {
        Matrix<T> output {matrix_a};
        T const temp = matrix_b(0,0);
        output = output * temp;
        return output;
    }

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


//Multiplies all elements of the matrix pair by a constant
template<typename T>
Matrix_pair<T> operator * (T const c, Matrix_pair<T> const& matrix_a) noexcept
{
    Matrix_pair<T> output{matrix_a};

    for (auto i = 0; i < output.first.max_rows(); ++i) {
        for (auto j = 0; j < output.first.max_cols(); ++j) {
            output.first(i,j) *= c;
        }
    }

    for (auto i = 0; i < output.second.max_rows(); ++i) {
        for (auto j = 0; j < output.second.max_cols(); ++j) {
            output.second(i,j) *= c;
        }
    }

    return output;
}


//Multiplies all elements of the matrix pair by a constant
template<typename T>
Matrix_pair<T> operator * (Matrix_pair<T> const& matrix_a, T const c) noexcept
{
    Matrix_pair<T> output{matrix_a};
    output = c * output;
    return output;
}


//Multiplies all elements of the matrix system by a constant
template<typename T>
Matrix_system<T> operator * (T const c, Matrix_system<T> const& matrix_a) noexcept
{
    Matrix_system<T> output{matrix_a};

    for (auto i = 0; i < output.first.max_rows(); ++i) {
        for (auto j = 0; j < output.first.max_cols(); ++j) {
            output.first(i,j) *= c;
        }
    }

    for (auto i = 0; i < output.second.max_rows(); ++i) {
        for (auto j = 0; j < output.second.max_cols(); ++j) {
            output.second(i,j) *= c;
        }
    }

    for (auto i = 0; i < output.third.max_rows(); ++i) {
        for (auto j = 0; j < output.third.max_cols(); ++j) {
            output.third(i,j) *= c;
        }
    }

    return output;
}


//Multiplies all elements of the matrix system by a constant
template<typename T>
Matrix_system<T> operator * (Matrix_system<T> const& matrix_a, T const c) noexcept
{
    Matrix_system<T> output{matrix_a};
    output = c * output;
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


//Performs the specified operation between a particular matrix_system row and a row slice. E.G. row_slice_op(A, S, 3, '+'); adds
//to row 3 of the first, second and third matrices the first, second and third slices from S
template<typename T>
void row_slice_op(Matrix_system<T>& m_sys, Matrix_system<T>& row_slice, sz_t row, char op) noexcept
{
    if (m_sys.first.max_cols()  != row_slice.first.size()  ||
        m_sys.second.max_cols() != row_slice.second.size() ||
        m_sys.third.max_cols()  != row_slice.third.size())
    {
        std::cout<<"Error in row_slice_op(): Matrix(ces) row size(s) is/are different than slice size(s)\n";
        return;
    }

    if (op != '+' && op != '-' && op != '*' && op != '/' && op != '=') {
        std::cout<<"Error in row_slice_op(): Unknown operation\n";
        return;
    }

    row_slice_op(m_sys.first, row_slice.first, row, op);
    row_slice_op(m_sys.second, row_slice.second, row, op);
    row_slice_op(m_sys.third, row_slice.third, row, op);

    return;
}


//Performs the specified operation between a particular Matrix_pair row and a Matrix_pair row slice. E.G.
//row_slice_op(A, S, 3, '+'); adds to row 3 of the first and second matrices, the first and second slices from S
template<typename T>
void row_slice_op(Matrix_pair<T>& m_sys, Matrix_pair<T>& row_slice, sz_t row, char op) noexcept
{
    if (m_sys.first.max_cols()  != row_slice.first.size()  ||
        m_sys.second.max_cols() != row_slice.second.size())
    {
        std::cout<<"Error in row_slice_op(): Matrix(ces) row size(s) is/are different than slice size(s)\n";
        return;
    }

    if (op != '+' && op != '-' && op != '*' && op != '/' && op != '=') {
        std::cout<<"Error in row_slice_op(): Unknown operation\n";
        return;
    }

    row_slice_op(m_sys.first, row_slice.first, row, op);
    row_slice_op(m_sys.second, row_slice.second, row, op);

}


//Performs the specified operation between a particular matrix column and a column slice. E.G. col_slice_op(A, S, 3, '+'); adds
//to column 3 of matrix A the slice S
template<typename T>
void col_slice_op(Matrix<T>& matrix, Matrix<T>& col_slice, sz_t col, char op) noexcept
{
    if (matrix.max_rows() != col_slice.size()) {
        std::cout<<"Error in row_slice_op(): Matrix row size is different than slice size\n";
        return;
    }

    switch(op)
    {
    case '+':
        {
            for (auto i = 0; i < matrix.max_rows(); ++i) {
                matrix(i, col) = matrix(i, col) + col_slice[i];
            }
            return;
        }

    case '-':
        {
            for (auto i = 0; i < matrix.max_rows(); ++i) {
                matrix(i, col) = matrix(i, col) - col_slice[i];
            }
            return;
        }
    case '*':
        {
            for (auto i = 0; i < matrix.max_rows(); ++i) {
                matrix(i, col) = matrix(i, col) * col_slice[i];
            }
            return;
        }

    case '/':
        {
            for (auto i = 0; i < matrix.max_rows(); ++i) {
                if (!(col_slice[i]*col_slice[i] > 0 )) {
                    std::cout<<"Error in col_slice_op(): Division by zero\n";
                    return;
                }

                matrix(i, col) = matrix(i, col) / col_slice[i];
            }
            return;
        }

    case '=':
        {
            for (auto i = 0; i < matrix.max_rows(); ++i) {
                matrix(i, col) = col_slice[i];
            }
            return;
        }

    default:
        {
            std::cout<<"Error in col_slice_op(): Unknown operation\n";
            return;
        }
    }

    return;
}


//Performs the specified operation between a particular matrix_system column and a column slice. E.G. col_slice_op(A, S, 3, '+');
//adds to column 3 of the first, second and third matrices the first, second and third slices from col_slice
template<typename T>
void col_slice_op(Matrix_system<T>& m_sys, Matrix_system<T>& col_slice, sz_t col, char op) noexcept
{
    if (m_sys.first.max_rows()  != col_slice.first.size()  ||
        m_sys.second.max_rows() != col_slice.second.size() ||
        m_sys.third.max_rows()  != col_slice.third.size())
    {
        std::cout<<"Error in col_slice_op(): Matrix(ces) row size(s) is/are different than slice size(s)\n";
        return;
    }

    if (op != '+' && op != '-' && op!= '*' && op != '/' && op != '=') {
        std::cout<<"Error in col_slice_op(): Unknown operation\n";
        return;
    }

    col_slice_op(m_sys.first, col_slice.first, col, op);
    col_slice_op(m_sys.second, col_slice.second, col, op);
    col_slice_op(m_sys.third, col_slice.third, col, op);

    return;
}


//Performs the specified operation between a particular Matrix_pair row and a Matrix_pair row slice. E.G.
//row_slice_op(A, S, 3, '+'); adds to row 3 of the first and second matrices, the first and second slices from S
template<typename T>
void col_slice_op(Matrix_pair<T>& m_sys, Matrix_pair<T>& col_slice, sz_t col, char op) noexcept
{
    if (m_sys.first.max_rows()  != col_slice.first.size()  ||
        m_sys.second.max_rows() != col_slice.second.size())
    {
        std::cout<<"Error in col_slice_op(): Matrix(ces) row size(s) is/are different than slice size(s)\n";
        return;
    }

    if (op != '+' && op != '-' && op!= '*' && op != '/' && op != '=') {
        std::cout<<"Error in col_slice_op(): Unknown operation\n";
        return;
    }

    col_slice_op(m_sys.first, col_slice.first, col, op);
    col_slice_op(m_sys.second, col_slice.second, col, op);

    return;
}


//inverts the target matrix. Note: Division by 0 errors imply that the matrix may not be invertible
template<typename T>
Matrix<T> invert(Matrix<T> const& A_) noexcept
{
    //Check if it is a square matrix
    if (A_.max_cols()!=A_.max_rows()) {
        std::cout<<"Error in invert(): Not a square matrix (m != n). Returning default constructed matrix\n";
        Matrix<T> output{};
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

#endif // MATRIX_H_INCLUDED
