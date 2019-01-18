#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

//Include list
#include <iostream>
#include <exception>
#include <string_view>
#include <optional>
#include <vector>
#include <thread>
#include <future>
#include <chrono>
#include <array>
#include <cmath>


using sz_t = std::size_t;
using sz_t = std::size_t;

//Number of threads available in this machine
static std::atomic<std::size_t> THREADS {std::thread::hardware_concurrency()};


namespace Maths
{

    template<typename T = sz_t>
    constexpr std::vector<std::array<T,2>> split_interval(T const st, T const ed)
    {
        auto const distance {ed-st};
        auto const step_size{distance / THREADS};
        std::vector<std::array<T, 2>> output;
        auto first {st};
        auto second {st};
    
        for (auto i = 0; i < THREADS; ++i)
        {
            if (i == THREADS-1)
            {
                second = ed;
                output.emplace_back(std::array<T, 2>{first, second});
                return output;
            }
        
            second = first + step_size;
            output.emplace_back(std::array<T, 2>{first, second});
            first = second;
        }
    
        return output;
    } //End of split_interval

    template<typename T>
    bool is_zero(T const c)
    {
        T cmp_val {static_cast<T>(1.0e-10)};
        if (std::abs(c) < cmp_val) return true;
        return false;
    }
    
    
    
    //Generic wrapper to create math errors
    struct Error final : public std::exception
    {
        explicit Error() : error_msg{std::string{"Error state not set"}} {}
        explicit Error(std::string_view error_msg_) : error_msg{error_msg_} {}
        
        char const* what() const noexcept
        {
            return error_msg.data();
        }
        
        Error& operator = (std::string_view message) noexcept
        {
            error_msg = message;
            return *this;
        }
        
        std::string error_msg;
    };
    
    

//Stores data as a flat vector, with overloaded call operator to be able to access the data as a 0 indexed matrix
//e.g. A(0,1) returns the element on the first row and second column. Recommended types: float, double. It also supports
//other types but factorizations will no longer work properly
    template<typename T>
    class Matrix
    {
    private:
        std::vector<T> data; //Container of matrix data
        sz_t m; //Number of rows
        sz_t n; //Number of columns

    public:
        //Constructors
        //Default constructor, single zero matrix
        Matrix() noexcept:
        m{1}, n{1}
        {
            data.resize(1, static_cast<T>(0.0));
        }


        //Construct matrix of zeros given dimensions
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


        //Constructor moving the potentially large data_
        template <typename inner_T>
        Matrix(sz_t m_, sz_t n_, inner_T&& data_) noexcept:
        m {m_}, n {n_}, data {m_ * n_ == data_.size() ? std::forward<inner_T>(data_) : 
            throw Error{"Error in Matrix constructor: Vector size does not correspond to matrix size\n"}} {}


        Matrix(Matrix<T> const& mat) = default;


        //Move assignment operator
        Matrix<T>& operator = (Matrix<T>&& matrix) noexcept = default;


        //Copy assignment operator
        Matrix<T>& operator = (Matrix<T> const& matrix) = default;
    

        //Overloaded operators
        //Access an element by inputting the row and column value, constant version
        T const& operator () (sz_t const i, sz_t const j) const
        {
            if (i < m && j < n) {
                return data[i + (m * j)];
            } else {
                if (j < m) {
                    std::string const e {"Error in operator(): row index "+std::to_string(i)+" exceeds m-1\n"};
                    throw Error{e};
                } else {
                    std::string const e {"Error in operator(): column index "+std::to_string(j)+" exceeds n-1\n"};
                    throw Error{e};
                }
            }
        }


        //Access an element by inputting the row and column value, non-const version
        T& operator () (sz_t const i, sz_t const j)
        {
            if (i < m && j < n) {
                return data[i + (m * j)];
            } else {
                if (j < m) {
                    std::string const e {"Error in operator(): row index "+std::to_string(i)+" exceeds m-1\n"};
                    throw Error{e};
                } else {
                    std::string const e {"Error in operator(): column index "+std::to_string(j)+" exceeds n-1\n"};
                    throw Error{e};
                }
            }
        }
        

        //Access the vector that contains the elements directly by index const version
        T const& operator [] (sz_t const i) const
        {
            if (i < data.size()) {
                return data[i];
            } else {
                std::string const e{"Error in operator[]: Index "+std::to_string(i)+" exceeds maximum size\n"};
                throw Error{e};
            }
        }


        //Access the vector that contains the elements directly by index non-const version
        T& operator [] (sz_t const i)
        {
            if (i < data.size()) {
                return data[i];
            } else {
                std::string const e{"Error in operator[]: Index "+std::to_string(i)+" exceeds maximum size\n"};
                throw Error{e};
            }
        }
        

        //Utility functions
        //Return a column slice or row slice. Note: This slice is a copy and not part of the matrix
        Matrix<T> slice (sz_t const k, char c) const
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
                        std::string const e {"Error in Matrix.slice(): Invalid column\n"};
                        throw Error{e};
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
                        std::string const e {"Error in Matrix.slice(): Invalid row\n"};
                        throw Error{e};
                    }
                }

            default:
                {
                    std::string e {"Error in Matrix.slice(): Unknown operation\n"};
                    throw Error{e};
                }
            }
        }


        //Return the size of the underlying vector container
        constexpr std::size_t size() const noexcept
        {
            return data.size();
        }


        //Change the underlying vector
        void set_data(std::vector<T>& data_) noexcept
        {
            std::swap(data_, data);
            return;
        }


        //Returns the underlying vector
        constexpr std::vector<T> const& get_data() const noexcept
        {
            return data;
        }
    
        //Makes a copy of the underlying vector
        std::vector<T> copy_vector() const noexcept
        {
            std::vector<T> out{data};
            return out;
        }


        //Returns the number of rows in the matrix
        constexpr sz_t max_rows() const noexcept
        {
            return m;
        }


        //Returns the number of columns in the matrix
        constexpr sz_t max_cols() const noexcept
        {
            return n;
        }


        //Changes the number of rows in the vector, either adding zeros at the end (increase) or truncating the matrix (decrease)
        void set_rows(sz_t m_, T init_val = 0)
        {
            if (m == m_) return;
            if (m_ == 0) {
                std::string const e {"Error in Matrix.set_rows(): Number of rows cannot be 0\n"};
                throw Error{e};
            }

            sz_t difference = m > m_ ? m-m_:m_-m;
            std::vector<T> new_vector;
            new_vector.resize(m_*n, init_val);
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
            std::swap(data, new_vector);
            return;
        } //End of set rows


        //Changes the number of columns in the vector, either adding zeros at the end (increase) or truncating the matrix (decrease)
        void set_cols(sz_t n_, T init_val = 0)
        {
            if (n == n_) return;
            if (n_ == 0) {
                std::string const e {"Error in Matrix.set_rows(): Number of columns cannot be 0\n"};
                throw Error{e};
            }

            std::vector<T> new_vector;
            new_vector.resize(m*n_, init_val);

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
            std::swap(data, new_vector);
            return;
        } //End of set cols

        friend void swap(Matrix& lhs, Matrix& rhs) noexcept
        {
            std::swap(lhs.data, rhs.data);
            std::swap(lhs.m, rhs.m);
            std::swap(lhs.n, rhs.n);
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

    }; //End of matrix


    //Helper structs
    //Auxiliary struct for returning pair of matrices
    template<typename T>
    struct Matrix_pair 
    {
        Matrix<T> first;
        Matrix<T> second;

        //Constructors
        //Default constructor
        Matrix_pair() noexcept:
        first {Matrix<T>()}, second {Matrix<T>()} {};

        //Copy constructor
        Matrix_pair(Matrix<T> const& first_, Matrix<T> const& second_) :
        first {first_}, second {second_} {}
        
        //Copy constructor2
        Matrix_pair(Matrix_pair<T> const& m_pair) = default;
        
        //Move constructor
        Matrix_pair(Matrix<T>&& first_, Matrix<T>&& second_) noexcept
        {
            swap(first_, first);
            swap(second_, second);
        }
        
        //Move constructor2
        Matrix_pair(Matrix_pair<T>&& m_pair) noexcept
        {
            swap(m_pair.first, first);
            swap(m_pair.second, second);
        }
        
        //Copy assignment
        Matrix_pair<T>& operator = (Matrix_pair<T> const& m_pair) = default;
        
        //Move assignment
        Matrix_pair<T>& operator = (Matrix_pair<T>&& m_pair) noexcept
        {
            swap(m_pair.first, first);
            swap(m_pair.second, second);
            return *this;
        }
        
        
        //Returns a column slice or row slice of the entire system. Note: This slice is a copy and not part of the matrix
        Matrix_pair<T> slice (sz_t const k, char c) const
        {

            switch(c)
            {
            case 'c':
                {
                    Matrix<T> A = first.slice(k, 'c');
                    Matrix<T> B = second.slice(k, 'c');
                    Matrix_pair<T> output {std::move(A), std::move(B)};
                    return output;
                }

            case 'r':
                {
                    Matrix<T> A = first.slice(k, 'r');
                    Matrix<T> B = second.slice(k, 'r');
                    Matrix_pair<T> output {std::move(A), std::move(B)};
                    return output;
                }

            default:
                {
                    std::string e{"Error in Matrix_pair.slice(): Unknown option\n"};
                    throw Error{e};
                }
            }
        }
    }; //Matrix_pair end


    //Auxiliary struct to represent system of matrices
    template<typename T>
    struct Matrix_system 
    {
        Matrix<T> first;
        Matrix<T> second;
        Matrix<T> third;


        //Constructors
        //Default constructor
        Matrix_system() noexcept:
        first {Matrix<T>()}, second {Matrix<T>()}, third {Matrix<T>{}} {};

        //Copy constructor
        Matrix_system(Matrix<T> const& first_, Matrix<T> const& second_, Matrix<T> const& third_):
        first {first_}, second {second_}, third {third_} {}
        
        //Copy constructor2
        Matrix_system(Matrix_system<T> const& m_pair) = default;
        
        //Move constructor
        Matrix_system(Matrix<T>&& first_, Matrix<T>&& second_, Matrix<T>&& third_) noexcept
        {
            swap(first_, first);
            swap(second_, second);
            swap(third_, third);
        }
        
        //Move constructor2
        Matrix_system(Matrix_system<T>&& m_system) noexcept
        {
            swap(m_system.first, first);
            swap(m_system.second, second);
            swap(m_system.third, third);
        }
        
        //Copy assignment
        Matrix_system<T>& operator = (Matrix_system<T> const& m_pair) = default;
        
        //Move assignment
        Matrix_system<T>& operator = (Matrix_system<T>&& m_system) noexcept
        {
            swap(m_system.first, first);
            swap(m_system.second, second);
            swap(m_system.third, third);
            return *this;
        }

        //Returns a column slice or row slice of the entire system. Note: This slice is a copy and not part of the matrix
        Matrix_system<T> slice (sz_t const k, char c) const
        {

            switch(c)
            {
            case 'c':
                {
                    Matrix<T> A = first.slice(k, 'c');
                    Matrix<T> B = second.slice(k, 'c');
                    Matrix<T> C = third.slice(k, 'c');
                    Matrix_system<T> output {std::move(A), std::move(B), std::move(C)};
                    return output;
                }

            case 'r':
                {
                    Matrix<T> A = first.slice(k, 'r');
                    Matrix<T> B = second.slice(k, 'r');
                    Matrix<T> C = third.slice(k, 'r');
                    Matrix_system<T> output {std::move(A), std::move(B), std::move(C)};
                    return output;
                }

            default:
                {
                    std::string const e {"Error in Matrix_system.slice(): Unknown option\n"};
                    throw Error{e};
                }

            }
        }
    }; //Matrix system


    //Other operators
    //Overloading << operator to help print the matrix or store it in a file.
    template<typename T>
    std::ostream& operator << (std::ostream& os, Matrix<T> const& matrix) noexcept
    {
        auto max_r {matrix.max_rows()};
        auto max_c {matrix.max_cols()};
        for (auto i = 0; i < max_r; ++i) {
            for (auto j = 0; j < max_c; ++j) {
                if (j != max_c-1) {
                    os<<matrix(i,j)<<' ';
                } else {
                os<<matrix(i,j)<<'\n';
                }
            }
        }

        return os;
    }


    //Overloaded << operator to print a matrix pair
    template<typename T>
    std::ostream& operator << (std::ostream& os, Matrix_pair<T> const& matrix) noexcept
    {
        std::cout<<"First matrix\n";
        os<<matrix.first<<'\n';

        std::cout<<"Second matrix\n";
        os<<matrix.second<<'\n';

        return os;
    }


    //Overloaded << operator to print a matrix system
    template<typename T>
    std::ostream& operator << (std::ostream& os, Matrix_system<T> const& matrix) noexcept
    {
        std::cout<<"First matrix\n";
        os<<matrix.first<<'\n';

        std::cout<<"Second matrix\n";
        os<<matrix.second<<'\n';

        std::cout<<"Third matrix\n";
        os<<matrix.third<<'\n';

        return os;
    }

    
    //Operations
    //Adds a constant to every element of the matrix,
    template<typename T>
    Matrix<T> operator + (Matrix<T> const& matrix_a, T const c) noexcept
    {
        auto max_r {matrix_a.max_rows()};
        auto max_c {matrix_a.max_cols()};
        Matrix<T> output {matrix_a};
        
        for (auto i = 0; i < max_r; ++i) {
            for (auto j = 0; j < max_c; ++j) {
                output(i,j) += c;
            }
        }
            
        return output;
    }


    //Adds a constant to every element of the matrix
    template<typename T>
    Matrix<T> operator + (T const c, Matrix<T> const& matrix_a) noexcept
    {
        Matrix<T> output {(matrix_a + c)};
        return output;
    }


    //Adds matrix_b to matrix_a
    template<typename T>
    Matrix<T> operator + (Matrix<T> const& matrix_a, Matrix<T> const& matrix_b)
    {
        Matrix<T> output {matrix_a};
        auto max_r {matrix_a.max_rows()};
        auto max_c {matrix_a.max_cols()};

        if (matrix_a.max_rows() != matrix_b.max_rows() || matrix_a.max_cols() != matrix_b.max_cols()) {
            std::string const e {"Error in Matrix operator +: Incompatible matrix dimensions\n"};
            throw Error{e};
        }

        for (auto i = 0; i < max_r; ++i) {
            for (auto j = 0; j < max_c; ++j) {
                output(i,j) = matrix_a(i,j) + matrix_b(i,j);
            }
        }

        return output;
    }


    //Substracts a constant to all elements of the matrix, the output is a new matrix
    template<typename T>
    Matrix<T> operator - (Matrix<T> const& matrix_a, T const c) noexcept
    {
        Matrix<T> output {matrix_a};
        auto max_r {matrix_a.max_rows()};
        auto max_c {matrix_a.max_cols()};
        
        for (auto i = 0; i < max_r; ++i) {
            for (auto j = 0; j < max_c; ++j) {
                output(i,j) -= c;
            }
        }

        return output;
    }


    //Substracts a constant to all elements of the matrix, the output is a new matrix
    template<typename T>
    Matrix<T> operator - (T const c, Matrix<T> const& matrix_a) noexcept
    {
        Matrix<T> output {matrix_a};
        
        auto max_r {matrix_a.max_rows()};
        auto max_c {matrix_a.max_cols()};
        
        for (auto i = 0; i < max_r; ++i) {
            for (auto j = 0; j < max_c; ++j) {
                output(i,j) = c - matrix_a(i,j);
            }
        }

        return output;
    }


    //Substracts matrix_b from matrix_a
    template<typename T>
    Matrix<T> operator - (Matrix<T> const& matrix_a, Matrix<T> const& matrix_b)
    {
        Matrix<T> output {matrix_a};
        auto max_r {matrix_a.max_rows()};
        auto max_c {matrix_a.max_cols()};
        
        if (matrix_a.max_rows() != matrix_b.max_rows() || matrix_a.max_cols() != matrix_b.max_cols()) {
            std::string const e {"Error in Matrix operator -: Incompatible matrix dimensions\n"};
            throw Error{e};
        }

        for (auto i = 0; i < max_r; ++i) {
            for (auto j = 0; j < max_c; ++j) {
                output(i, j) = matrix_a(i,j) - matrix_b(i,j);
            }
        }

        return output;
    }


    //Multiplies all elements of the matrix by a constant
    template<typename T>
    Matrix<T> operator * (Matrix<T> const& matrix_a, T const c) noexcept
    {
        Matrix<T> output {matrix_a};
        auto max_r {matrix_a.max_rows()};
        auto max_c {matrix_a.max_cols()};
        for (auto i = 0; i < max_r; ++i) {
            for (auto j = 0; j < max_c; ++j) {
                output(i,j) *= c;
            }
        }

        return output;
    }


    //Multiplies all elements of the matrix by a constant
    template<typename T>
    Matrix<T> operator * (T const c, Matrix<T> const& matrix_a) noexcept
    {
        Matrix<T> output {matrix_a};
        output = output * c; 

        return output;
    }
    
    
    //Multiplies two matrices of compatible dimensions.
    template<typename T>
    Matrix<T> operator * (Matrix<T> const& matrix_a, Matrix<T> const& matrix_b)
    {
        if (matrix_a.max_cols() != matrix_b.max_rows()) {
            std::string const e {"Error in Matrix operator *: Incompatible matrix dimensions\n"};
            throw Error{e};
        }

        sz_t m_ = matrix_a.max_rows();
        sz_t n_a = matrix_a.max_cols();
        sz_t n_ = matrix_b.max_cols();
        Matrix<T> output(m_, n_);

        for (auto i = 0; i < m_; ++i) {
            for (auto j = 0; j < n_; ++j) {
                for (auto k = 0; k < n_a; ++k) {
                    output(i, j) += matrix_a(i, k) * matrix_b(k, j);
                }
            }
        }

        return output;
    }


    //Divides all elements of the matrix by a constant
    template<typename T>
    Matrix<T> operator / (Matrix<T> const& matrix_a, T const c)
    {
        if (is_zero<T>(c))
        {
            std::string const e {"Error in Matrix operator /: Division by zero encountered\n"};
            throw Error{e};
        }
        
        Matrix<T> output {matrix_a};
        auto max_r {matrix_a.max_rows()};
        auto max_c {matrix_a.max_cols()};

        for (auto i = 0; i < max_r; ++i) {
            for (auto j = 0; j < max_c; ++j) {
                output(i,j) /= c;
            }
        }

        return output;
    }
    
    
    //Matrix pair operations, Note they don't return a new instance of matrix pair
    //Multiplies all elements of the matrix pair by a constant
    template<typename T>
    void operator * (T const c, Matrix_pair<T>& matrix_a) noexcept
    {
        matrix_a.first = c * matrix_a.first; 
        matrix_a.second = c * matrix_a.second; 
        return;
    }


    //Multiplies all elements of the matrix pair by a constant
    template<typename T>
    void operator * (Matrix_pair<T>& matrix_a, T const c) noexcept
    {
        c * matrix_a;
        return;
    }


    //Multiplies all elements of the matrix system by a constant
    template<typename T>
    void operator * (T const c, Matrix_system<T>& matrix_a) noexcept
    {

        matrix_a.first = matrix_a.first * c;
        matrix_a.second = matrix_a.second * c;
        matrix_a.third = matrix_a.third * c;

        return;
    }


    //Multiplies all elements of the matrix system by a constant
    template<typename T>
    void operator * (Matrix_system<T>& matrix_a, T const c) noexcept
    {
        c * matrix_a;
        return;
    }


    //Utility functions
    //Sets all available threads to work on a function.
    template<typename f_ptr, typename...ARGS>
    void launch_task (f_ptr const F, 
                      sz_t const starts, 
                      sz_t const ends,
                      int const delay_ms, 
                      ARGS...args) noexcept
    {
        std::vector<std::array<sz_t, 2>> intervals {split_interval(starts, ends)};
        
        std::vector<std::thread> threads;
        for (auto i = 0; i < THREADS; ++i) {
            std::thread t {F, intervals[i][0], intervals[i][1], args...};
            threads.emplace_back<std::thread>(std::move(t));
            std::this_thread::sleep_for(std::chrono::milliseconds(delay_ms));
        }

        for (auto i = 0; i < threads.size(); ++i) {
            threads[i].join();
        }

        return;
    }


    //Switches the rows with the columns of the matrix
    template<typename T>
    Matrix<T> transpose(Matrix<T> const& matrix_a, sz_t const limit = 500000) noexcept
    {
        auto max_r {matrix_a.max_rows()};
        auto max_c {matrix_a.max_cols()};
        auto out_r {max_c};
        auto out_c {max_r};
        
        Matrix<T> output {out_r, out_c};
        
        auto inner_launcher = [&] (sz_t const st, sz_t const ed, bool col_split = false) 
        {
            //Check whether to split by rows or columns
            if (col_split)
            {
                for (auto i = 0; i < out_r; ++i)
                {
                    for (auto j = st; j < ed; ++j)
                    {
                        output(i,j) = matrix_a(j, i);
                    }
                }
            }
            else
            {
                for (auto i = st; i < ed; ++i)
                {
                    for (auto j = 0; j < out_c; ++j)
                    {
                        output(i,j) = matrix_a(j, i);
                    }
                }
            }
            return;
        };
        
        if (matrix_a.size() > limit && limit > THREADS)
        {
            //Split by columns if they are larger, else split by rows
            if (out_r < out_c)
            {
                launch_task(inner_launcher, 0, out_c, 0, true);
            } 
            else
            {
                launch_task(inner_launcher, 0, out_r, 0, false);
            }
            
        } 
        else 
        {
            inner_launcher(0, out_r);
        }
        
        return output;
    } //End of transpose


    //Creates a matrix with 0 in all the elements except those in the diagonal wich are equal to 1
    template<typename T>
    Matrix<T> make_identity(sz_t const m, sz_t const n) noexcept
    {
        Matrix<T> output {m, n};
        for (auto i = 0; i < m; ++i)
        {
            output(i,i) = static_cast<T>(1.0);
        }
    
        return output;
    }


    //Performs the specified operation between a particular matrix row and a row slice. E.G. row_slice_op(A, S, 3, '+'); adds
    //to row 3 of matrix A the slice S
    template<typename T>
    Matrix<T> row_slice_op(Matrix<T> const& matrix, Matrix<T> const& row_slice, sz_t const row, char const op)
    {
        if (matrix.max_cols() != row_slice.size()) 
        {
            std::string const e {"Error in row_slice_op(): Matrix row size is different than slice size\n"};
            throw Error{e};
        }
        
        Matrix<T> output {matrix};
        auto max_c {matrix.max_cols()};
        try
        {
            switch(op)
            {
            case '+':
                {
                    for (auto j = 0; j < max_c; ++j) {
                        output(row, j) = output(row, j) + row_slice[j];
                    }
                
                    return output;
                }

            case '-':
                {
                    for (auto j = 0; j < max_c; ++j) {
                        output(row, j) = output(row, j) - row_slice[j];
                    }
                
                    return output;
                }
            case '*':
                {
                    for (auto j = 0; j < max_c; ++j) {
                    output(row, j) = output(row, j) * row_slice[j];
                    }
                
                    return output;
                }

            case '/':
                {
                    for (auto j = 0; j < max_c; ++j) {
                        if (is_zero<T>(row_slice[j]))
                        {
                            std::string const e {"Error in row_slice_op(): Division by zero\n"};
                            throw Error{e};
                        }
                    
                        output(row, j) = output(row, j) / row_slice[j];
                    }
                
                    return output;
                }

            case '=':
                {
                    for (auto j = 0; j < max_c; ++j) {
                        output(row, j) = row_slice[j];
                    }
                
                    return output;
                }

            default:
                {
                    std::string const e {"Error in row_slice_op(): Unknown operation\n"};
                    throw Error{e};
                }
            }
        }
        catch(Error& e)
        {
            std::string const prev_e {e.what()};
            std::string const curr_e {"Error in row_slice_op() -> "+prev_e};
            throw Error{curr_e};
        }
    } //end row_slice_op


    //Performs the specified operation between a particular Matrix_pair row and a Matrix_pair row slice. E.G.
    //row_slice_op(A, S, 3, '+'); adds to row 3 of the first and second matrices, the first and second slices from S
    //Rethrows exception from single matrix row_slice_op
    template<typename T>
    void row_slice_op(Matrix_pair<T>& m_sys, Matrix_pair<T> const& row_slice, sz_t const row, char const op)
    {
        try
        {
            m_sys.first  = std::move(row_slice_op(m_sys.first, row_slice.first, row, op));
            m_sys.second = std::move(row_slice_op(m_sys.second, row_slice.second, row, op));
        }
        catch (std::runtime_error& e)
        {
            std::string const caught_e {e.what()};
            std::string const this_e {"Error in row_slice_op (Matrix_pair)-> "+caught_e};
            throw Error{this_e};
        }
        return;
    }
    

    //Performs the specified operation between a particular matrix_system row and a row slice. E.G. row_slice_op(A, S, 3, '+'); adds
    //to row 3 of the first, second and third matrices the first, second and third slices from S
    //Rethrows exception from single matrix row_slice_op
    template<typename T>
    void row_slice_op(Matrix_system<T>& m_sys, Matrix_system<T> const& row_slice, sz_t const row, char const op)
    {
        try
        {
            m_sys.first  = std::move(row_slice_op(m_sys.first , row_slice.first , row, op));
            m_sys.second = std::move(row_slice_op(m_sys.second, row_slice.second, row, op));
            m_sys.third  = std::move(row_slice_op(m_sys.third , row_slice.third , row, op));
        }
        catch (std::runtime_error& e)
        {
            std::string const caught_e {e.what()};
            std::string const this_e {"Error in row_slice_op (Matrix_system)-> "+caught_e};
            throw Error{this_e};
        }
        return;
    } //end of row_slice_op 


    //Performs the specified operation between a particular matrix column and a column slice. E.G. col_slice_op(A, S, 3, '+'); adds
    //to column 3 of matrix A the slice S
    //Allows exception handling
    template<typename T>
    Matrix<T> col_slice_op(Matrix<T> const& matrix, Matrix<T> const& col_slice, sz_t const col, char const op)
    {
        if (matrix.max_rows() != col_slice.size()) {
            std::string const e {"Error in col_slice_op(): Matrix row size is different than slice size\n"};
            throw Error{e};
        }
        
        Matrix<T> output {matrix};
        auto max_r {matrix.max_rows()};
        
        try
        {
            switch(op)
            {
            case '+':
                {
                    for (auto i = 0; i < max_r; ++i) {
                        output(i, col) = output(i, col) + col_slice[i];
                    }
                
                    return output;
                }

            case '-':
                {
                    for (auto i = 0; i < max_r; ++i) {
                        output(i, col) = output(i, col) - col_slice[i];
                    }
                
                    return output;
                }
            case '*':
                {
                    for (auto i = 0; i < max_r; ++i) {
                        output(i, col) = output(i, col) * col_slice[i];
                    }
                
                    return output;
                }

            case '/':
                {
                    for (auto i = 0; i < max_r; ++i) {
                        if (is_zero<T>(col_slice[i])) {
                            std::string const e{"Error in col_slice_op(): Division by zero\n"};
                            throw Error{e};
                        }

                        output(i, col) = output(i, col) / col_slice[i];
                    }
                
                    return output;
                }

            case '=':
                {
                    for (auto i = 0; i < max_r; ++i) {
                        output(i, col) = col_slice[i];
                    }
                
                    return output;
                }

            default:
                {
                    std::string const e {"Error in col_slice_op(): Unknown operation\n"};
                    throw Error{e};
                }
            }
        }
        catch (Error& e)
        {
            std::string const prev_e {e.what()};
            std::string const curr_e {"Error in col_slice_op() -> "+prev_e};
            throw Error{curr_e};
        }
    } //end of col_slice_op


    //Performs the specified operation between a particular Matrix_pair row and a Matrix_pair row slice. E.G.
    //row_slice_op(A, S, 3, '+'); adds to row 3 of the first and second matrices, the first and second slices from S
    template<typename T>
    void col_slice_op(Matrix_pair<T>& m_sys, Matrix_pair<T> const& col_slice, sz_t const col, char const op)
    {
        try
        {
            m_sys.first  = std::move(col_slice_op(m_sys.first , col_slice.first , col, op));
            m_sys.second = std::move(col_slice_op(m_sys.second, col_slice.second, col, op));
        }
        catch(std::runtime_error& e)
        {
            std::string const prev_e {e.what()};
            std::string const curr_e {"Error in col_slice_op (Matrix_pair) ->: "+prev_e};
            throw  Error{curr_e};
        }
        
        return;
    }
    
    
    //Performs the specified operation between a particular matrix_system column and a column slice. E.G. col_slice_op(A, S, 3, '+');
    //adds to column 3 of the first, second and third matrices the first, second and third slices from col_slice
    //Allows exception handling
    template<typename T>
    void col_slice_op(Matrix_system<T>& m_sys, Matrix_system<T> const& col_slice, sz_t const col, char const op)
    {
        try
        {
            m_sys.first  = std::move(col_slice_op(m_sys.first , col_slice.first , col, op));
            m_sys.second = std::move(col_slice_op(m_sys.second, col_slice.second, col, op));
            m_sys.third  = std::move(col_slice_op(m_sys.thrid , col_slice.third , col, op));
        }
        catch(std::runtime_error& e)
        {
            std::string const prev_e {e.what()};
            std::string const curr_e {"Error in col_slice_op (Matrix_system) -> "+prev_e};
            throw Error {curr_e};
        }
        
        return;
    } //End of col_slice_op


} //End namespace Maths

#endif