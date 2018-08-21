#ifndef FLAT_VEC_H_INCLUDED
#define FLAT_VEC_H_INCLUDED
#include <vector>
#include <iostream>


using sz_t = std::size_t;


template<typename T>
class Flat_vec
{
private:
    sz_t I;
    sz_t J;
    sz_t K;
    std::vector<T> vec;

public:
    //Constructors
    //Default constructor
    Flat_vec() noexcept:
        I{1}, J{1}, K{1}
    {
        vec.resize(I*J*K);
    }


    explicit Flat_vec(sz_t const I_) noexcept
    {
        auto I_val = I_==0 ? 1:I_;
        I = I_val;
        J = 1;
        K = 1;
        vec.resize(I*J*K);
    }


    Flat_vec(sz_t const I_, sz_t const J_) noexcept
    {
        auto I_val = I_==0 ? 1:I_;
        auto J_val = J_==0 ? 1:J_;
        I = I_val;
        J = J_val;
        K = 1;
        vec.resize(I*J*K);
    }


    Flat_vec(sz_t const I_, sz_t const J_, sz_t const K_) noexcept
    {
        auto I_val = I_==0 ? 1:I_;
        auto J_val = J_==0 ? 1:J_;
        auto K_val = K_==0 ? 1:K_;
        I = I_val;
        J = J_val;
        K = K_val;
        vec.resize(I*J*K);
    }


    //Copy constructor
    Flat_vec(Flat_vec<T> const& a) noexcept
    {
        vec = a.get_vec();
        I = a.max_I();
        J = a.max_J();
        K = a.max_K();
    }


    //Move constructor
    Flat_vec(Flat_vec<T>&& a) = default;


    //Overloaded operators and non-member functions
    //Copy assignment
    Flat_vec<T>& operator = (Flat_vec<T> const& a) noexcept
    {
        vec = a.get_vec();
        I = a.max_I();
        J = a.max_J();
        K = a.max_K();
        return *this;
    }


    //Move assignment
    Flat_vec<T>& operator = (Flat_vec<T>&& a) = default;


    //Access operator, const version
    T const& operator [] (sz_t const i) const noexcept
    {
        if (i < vec.size()) {
            return vec[i];
        } else {
            std::cout<<"Error in operator []: index "<<i<<" exceeds max index "<<
            (vec.size()-1)<<", returning first element...\n";
            return vec[0];
        }
    }


    //Access operator
    T& operator [] (sz_t const i) noexcept
    {
        if (i < vec.size()) {
            return vec[i];
        } else {
            std::cout<<"Error in operator []: index "<<i<<" exceeds max index "<<
            (vec.size()-1)<<", returning first element...\n";
            return vec[0];
        }
    }


    //Call operator for accessing the elements in an x,y,z manner, const version
    T const& operator () (sz_t i, sz_t j = 0, sz_t k = 0) const noexcept
    {
        if ((i < I) && (j < J) && (k < K)) {
            return vec[i + (I * j) + (I * J * k)];
        } else {
            sz_t mI = max_I()> 0 ? max_I() - 1.0 : 0;
            sz_t mJ = max_J()> 0 ? max_J() - 1.0 : 0;
            sz_t mK = max_K()> 0 ? max_K() - 1.0 : 0;
            bool eI = i > mI;
            bool eJ = j > mJ;
            bool eK = k > mK;
            std::cout<<"Error in operator (): "<<(eI ? "i ":"")<<(eJ ? "j ":"")<<(eK ? "k ":"")<<
            "index/indexes exceed maximum values (i = "<<(mI)<<", j = "<<(mJ)<<", k = "<<(mK)<<
            ").\nReturning first element...\n";
            return vec[0];
        }
    }


    //Call operator for accessing the elements in an x,y,z manner
    T& operator () (sz_t i, sz_t j = 0, sz_t k = 0) noexcept
    {
        if ((i < I) && (j < J) && (k < K)) {
            return vec[i + (I * j) + (I * J * k)];
        } else {
            sz_t mI = max_I()> 0 ? max_I() - 1.0 : 0;
            sz_t mJ = max_J()> 0 ? max_J() - 1.0 : 0;
            sz_t mK = max_K()> 0 ? max_K() - 1.0 : 0;
            bool eI = i > mI;
            bool eJ = j > mJ;
            bool eK = k > mK;
            std::cout<<"Error in operator (): "<<(eI ? "i ":"")<<(eJ ? "j ":"")<<(eK ? "k ":"")<<
            "index/indexes exceed maximum values (i = "<<(mI)<<", j = "<<(mJ)<<", k = "<<(mK)<<
            ").\nReturning first element...\n";
            return vec[0];
        }
    }


    //Returns the size of the underlying vector
    unsigned long long size() const noexcept
    {
        return vec.size();
    }


    //Returns an iterator to the beginning of the vector
    auto& begin() noexcept
    {
        return vec.begin();
    }


    //Returns a constant iterator to the beginning of the vector
    auto& cbegin() const noexcept {
        return vec.cbegin();
    }


    //Returns an iterator to the end of the vector
    auto& end() noexcept {
        return vec.end();
    }


    //Returns a constant iterator to the end of the vector
    auto& cend() const noexcept {
        return vec.cend();
    }


    //Member functions
    //Returns the size of the I dimension
    sz_t max_I() const noexcept
    {
        return I;
    }


    void set_max_I(sz_t const new_I) noexcept
    {
        if (new_I == 0)
        {
            std::cout<<"Error in set_max_I(): Size cannot be less than 1\n";
            return;
        }

        I = new_I;
        vec.resize(I*J*K);
    }


    //Returns the size of the J dimension
    sz_t max_J() const noexcept
    {
        return J;
    }


    void set_max_J(sz_t const new_J) noexcept
    {
        if (new_J == 0)
        {
            std::cout<<"Error in set_max_J(): Size cannot be less than 1\n";
            return;
        }

        J = new_J;
        vec.resize(I*J*K);
    }


    //Returns the size of the K dimension
    sz_t max_K() const noexcept
    {
        return K;
    }


    //Changes the maximum number of rows
    void set_max_K(sz_t new_K) noexcept
    {
        if (new_K == 0)
        {
            std::cout<<"Error in set_max_K(): Size cannot be less than 1\n";
            return;
        }

        K = new_K;
        vec.resize(I*J*K);
    }


    //Changes the underlying vector
    void set_vec (std::vector<T> new_vec) noexcept
    {
        vec = new_vec;
        return;
    }


    //Returns a copy of the underlying vector
    std::vector<T> get_vec() const noexcept
    {
        return vec;
    }
};


//Overloaded operators
//Add two flat vecs
template<typename T>
Flat_vec<T> operator + (Flat_vec<T> const& a, Flat_vec<T> const& b) noexcept
{
        if (a.size() != b.size()) {
            std::cout<<"Error in operator +: a.size() = "<<a.size()<<" != b.size() = "<<b.size()<<
            ". Returning default constructed Flat_vec\n";
            Flat_vec<T> output{};
            return output;
        }

        Flat_vec<T> output{a};

        for (unsigned long long i = 0; i < a.size(); ++i) {
            output[i] = a[i] + b[i];
        }

        return output;
}


//Add a constant to Flat_Vec
template<typename T>
Flat_vec<T> operator + (Flat_vec<T> const& a, T const b)
{
    Flat_vec<T> output {a};
    for (auto i = 0; i < a.size(); ++i) {
        output[i] += b;
    }
    return output;
}


//Add a Flat_vec to a constant
template<typename T>
Flat_vec<T> operator + (T const b, Flat_vec<T> const& a)
{
    Flat_vec<T> output {a};
    for (auto i = 0; i < a.size(); ++i) {
        output[i] += b;
    }
    return output;
}


//Subtract two flat vecs
template<typename T>
Flat_vec<T> operator - (Flat_vec<T> const& a, Flat_vec<T> const& b) noexcept
{
        if (a.size() != b.size()) {
            std::cout<<"Error in operator -: a.size() = "<<a.size()<<" != b.size() = "<<b.size()<<'\n';
            Flat_vec<T> output{};
            return output;
        }

        Flat_vec<T> output{a};

        for (unsigned long long i = 0; i < a.size(); ++i) {
            output[i] = a[i] - b[i];
        }

        return output;
}


//Subtract a constant to Flat_Vec
template<typename T>
Flat_vec<T> operator - (Flat_vec<T> const& a, T const b)
{
    Flat_vec<T> output {a};
    for (auto i = 0; i < a.size(); ++i) {
        output[i] -= b;
    }
    return output;
}


//Subtract a Flat_vec to a constant
template<typename T>
Flat_vec<T> operator - (T const b, Flat_vec<T> const& a)
{
    Flat_vec<T> output {a};
    for (auto i = 0; i < a.size(); ++i) {
        output[i] += b;
    }
    return output;
}


//Multiply (element by element) two flat vecs
template<typename T>
Flat_vec<T> operator * (Flat_vec<T> const& a, Flat_vec<T> const& b) noexcept
{
        if (a.size() != b.size()) {
            std::cout<<"Error in operator /: a.size() = "<<a.size()<<" != b.size() = "<<b.size()<<
            ". Returning default constructed Flat_vec\n";
            Flat_vec<T> output{};
            return output;
        }

        Flat_vec<T> output{a};

        for (unsigned long long i = 0; i < a.size(); ++i) {
            output[i] = a[i] * b[i];
        }

        return output;
}


//Multiply a constant to Flat_Vec
template<typename T>
Flat_vec<T> operator * (Flat_vec<T> const& a, T const b)
{
    Flat_vec<T> output {a};
    for (auto i = 0; i < a.size(); ++i) {
        output[i] *= b;
    }
    return output;
}


//Multiply a Flat_vec to a constant
template<typename T>
Flat_vec<T> operator * (T const b, Flat_vec<T> const& a)
{
    Flat_vec<T> output {a};
    for (auto i = 0; i < a.size(); ++i) {
        output[i] *= b;
    }
    return output;
}


//Divide (element by element) two flat vecs
template<typename T>
Flat_vec<T> operator / (Flat_vec<T> const& a, Flat_vec<T> const& b) noexcept
{
        if (a.size() != b.size()) {
            std::cout<<"Error in operator /: a.size() = "<<a.size()<<" != b.size() = "<<b.size()<<
            ". Returning default constructed Flat_vec\n";
            Flat_vec<T> output{};
            return output;
        }

        Flat_vec<T> output{a};
        bool warn_flag = false;
        sz_t index = 0;

        for (unsigned long long i = 0; i < a.size(); ++i) {
            if (b[i] != 0) {
                output[i] = a[i] / b[i];
            } else {
                index = i;
                warn_flag = true;
                break;
            }

        }

        if (warn_flag == true) {
        std::cout<<"Warning in operator /: Division by zero encountered at index "<<index<<". Stopping operation\n";
        }

        return output;
}


//Divide a Flat_vec by a constant.
template<typename T>
Flat_vec<T> operator / (Flat_vec<T> const& a, T const b)
{
    if (b == 0) {
        std::cout<<"Warning in operator /: Division by zero encountered. Returning default constructed Flat_vec\n";
        Flat_vec<T> output{};
        return output;
    }

    Flat_vec<T> output {a};
    for (auto i = 0; i < a.size(); ++i) {
        output[i] /= b;
    }

    return output;
}


//Divide a Flat_vec by a constant
template<typename T>
Flat_vec<T> operator / (T const b, Flat_vec<T> const& a)
{
    if (b == 0) {
        std::cout<<"Warning in operator /: Division by zero encountered. Returning default constructed Flat_vec\n";
        Flat_vec<T> output{};
        return output;
    }

    Flat_vec<T> output {a};
    for (auto i = 0; i < a.size(); ++i) {
        output[i] /= b;
    }

    return output;
}
#endif // FLAT_VEC_H_INCLUDED
