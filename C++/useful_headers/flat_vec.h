#ifndef FLAT_VEC_H_INCLUDED
#define FLAT_VEC_H_INCLUDED
#include <vector>
#include <iostream>

using ullong = unsigned long long;


template<typename T>
class Flat_vec
{
public:
    //Constructors
    Flat_vec() noexcept {
        vec.resize(I*J*K);
    }

    Flat_vec(ullong const I) noexcept : I{I} {
        auto I_val = I==0 ? 1:I;
        vec.resize(I_val*J*K);
    }

    Flat_vec(ullong const I, ullong const J) noexcept : I{I}, J{J} {
        auto I_val = I==0 ? 1:I;
        auto J_val = J==0 ? 1:J;
        vec.resize(I_val*J_val*K);
    }

    Flat_vec(ullong const I, ullong const J, ullong const K) noexcept : I{I}, J{J}, K{K} {
        auto I_val = I==0 ? 1:I;
        auto J_val = J==0 ? 1:J;
        auto K_val = K==0 ? 1:K;
        vec.resize(I_val*J_val*K_val);
    }

    //Move constructor (Note does not delete the previous vector)
    Flat_vec(Flat_vec<T>&& a) noexcept : I{a.max_I()}, J {a.max_J()}, K {a.max_K()}, vec {a.return_vec()} {};

    //Copy constructor
    Flat_vec(Flat_vec<T>& a) noexcept : I{a.max_I()}, J {a.max_J()}, K {a.max_K()}, vec {a.return_vec()} {};


    //Overloaded operators and non-member functions
    //Move assignment

    Flat_vec<T>& operator = (Flat_vec<T>&& a) noexcept {
        set_max_I(a.max_I());
        set_max_J(a.max_J());
        set_max_K(a.max_K());
        set_vec(a.return_vec());
        return *this;
    }

    //Copy assignment
    Flat_vec<T>& operator = (Flat_vec<T>& a) noexcept {
        set_max_I(a.max_I());
        set_max_J(a.max_J());
        set_max_K(a.max_K());
        set_vec(a.return_vec());
        return *this;
    }


    T& operator [] (ullong const i) noexcept {
        if (i < vec.size()) {
            return vec[i];
        } else {
            std::cout<<"Error in operator [] : index "<<i<<" exceeds max index "<<
            (vec.size()-1)<<", returning first element...\n";
            return vec[0];
        }
    }

    T const & operator [] (ullong const i) const noexcept {
        if (i < vec.size()) {
            return vec[i];
        } else {
            std::cout<<"Error in operator [] : index "<<i<<" exceeds max index "<<
            (vec.size()-1)<<", returning first element...\n";
            return vec[0];
        }
    }

    T& operator () (ullong i, ullong j = 0, ullong k = 0) noexcept {
        if ((i < I) && (j < J) && (k < K)) {
            return vec[i + (I * j) + (I * J * k)];
        } else {
            double mI = static_cast<double> (max_I()) - 1.0;
            double mJ = static_cast<double> (max_J()) - 1.0;
            double mK = static_cast<double> (max_K()) - 1.0;
            bool eI = i > mI;
            bool eJ = j > mJ;
            bool eK = k > mK;
            std::cout<<"Error in operator () : "<<(eI ? "i ":"")<<(eJ ? "j ":"")<<(eK ? "k ":"")<<
            "index/indexes exceed maximum values (i = "<<(mI)<<", j = "<<(mJ)<<", k = "<<(mK)<<
            ").\nReturning first element...\n";
            return vec[0];
        }
    }

    T const & operator () (ullong i, ullong j = 0, ullong k = 0) const noexcept {
        if ((I == 0 || i < I) && (J ==0 || j < J) && (K == 0 || k < K)) {
            return vec[i + (I * j) + (I * J * k)];
        } else {
            double mI = static_cast<double> (max_I()) - 1.0;
            double mJ = static_cast<double> (max_J()) - 1.0;
            double mK = static_cast<double> (max_K()) - 1.0;
            bool eI = i > mI;
            bool eJ = j > mJ;
            bool eK = k > mK;
            std::cout<<"Error in operator () : "<<(eI ? "i ":"")<<(eJ ? "j ":"")<<(eK ? "k ":"")<<
            "index/indexes exceed maximum values (i = "<<(mI)<<", j = "<<(mJ)<<", k = "<<(mK)<<
            ").\nReturning first element...\n";
            return vec[0];
        }
    }

    void operator + (Flat_vec<T> const & b) noexcept {
        if (vec.size() != b.size()) {
            std::cout<<"Error in operator + : vec.size() = "<<vec.size()<<" != b.size() = "<<b.size()<<'\n';
            return;
        }

        for (unsigned long long i = 0; i < vec.size(); ++i) {
            vec[i] = vec[i] + b[i];
        }

        return;
    }

    void operator + (std::vector<T> const & b) noexcept {
        if (vec.size() != b.size()) {
            std::cout<<"Error in operator + : vec.size() = "<<vec.size()<<" != b.size() = "<<b.size()<<'\n';
            return;
        }

        for (unsigned long long i = 0; i < vec.size(); ++i) {
            vec[i] = vec[i] + b[i];
        }

        return;
    }

    void operator + (T const b) noexcept {
        for (unsigned long long i = 0; i < vec.size(); ++i) {
            vec[i] = vec[i] + b;
        }

        return;
    }

    void operator - (Flat_vec<T> const & b) noexcept {
        if (vec.size() != b.size()) {
            std::cout<<"Error in operator - : vec.size() = "<<vec.size()<<" != b.size() = "<<b.size()<<'\n';
            return;
        }

        for (unsigned long long i = 0; i < vec.size(); ++i) {
            vec[i] = vec[i] - b[i];
        }

        return;
    }

    void operator - (std::vector<T> const & b) noexcept {
        if (vec.size() != b.size()) {
            std::cout<<"Error in operator - : vec.size() = "<<vec.size()<<" != b.size() = "<<b.size()<<'\n';
            return;
        }

        for (unsigned long long i = 0; i < vec.size(); ++i) {
            vec[i] = vec[i] - b[i];
        }

        return;
    }

    void operator - (T const b) noexcept {
        for (unsigned long long i = 0; i < vec.size(); ++i) {
            vec[i] = vec[i] - b;
        }

        return;
    }

    void operator * (Flat_vec<T> const & b) noexcept {
        if (vec.size() != b.size()) {
            std::cout<<"Error in operator * : vec.size() = "<<vec.size()<<" != b.size() = "<<b.size()<<'\n';
            return;
        }

        for (unsigned long long i = 0; i < vec.size(); ++i) {
            vec[i] = vec[i] * b[i];
        }

        return;
    }

    void operator * (std::vector<T> const & b) noexcept {
        if (vec.size() != b.size()) {
            std::cout<<"Error in operator * : vec.size() = "<<vec.size()<<" != b.size() = "<<b.size()<<'\n';
            return;
        }

        for (unsigned long long i = 0; i < vec.size(); ++i) {
            vec[i] = vec[i] * b[i];
        }

        return;
    }

    void operator * (T const val) noexcept {

        for (unsigned long long i = 0; i < vec.size(); ++i) {
            vec[i] = vec[i] * val;
        }

        return;
    }

    void operator / (Flat_vec<T> const & b) noexcept {
        if (vec.size() != b.size()) {
            std::cout<<"Error in operator / : vec.size() = "<<vec.size()<<" != b.size() = "<<b.size()<<'\n';
            return;
        }

        for (unsigned long long i = 0; i < size(); ++i) {
            vec[i] = vec[i] / b[i];
        }

        return;
    }

    void operator / (std::vector<T> const & b) noexcept {
        if (vec.size() != b.size()) {
            std::cout<<"Error in operator / : vec.size() = "<<vec.size()<<" != b.size() = "<<b.size()<<'\n';
            return;
        }

        for (unsigned long long i = 0; i < size(); ++i) {
            vec[i] = vec[i] / b[i];
        }

        return;
    }

    void operator / (T const val) noexcept {
        for (unsigned long long i = 0; i < size(); ++i) {
            vec[i] = vec[i] / val;
        }

        return;
    }

    unsigned long long size() const noexcept {
        return vec.size();
    }

    auto& begin() noexcept {
        return vec.begin();
    }

    auto& begin() const noexcept {
        return vec.cbegin();
    }

    auto& end() noexcept {
        return vec.end();
    }

    auto& end() const noexcept {
        return vec.cend();
    }


    //Member functions
    constexpr ullong max_I() const noexcept {
        return I;
    }

    constexpr void set_max_I(ullong new_I) noexcept {
        I = new_I;
        return;
    }

    constexpr ullong max_J() const noexcept {
        return J;
    }

    constexpr void set_max_J(ullong new_J) noexcept {
        J = new_J;
        return;
    }

    constexpr ullong max_K() const noexcept {
        return K;
    }

    constexpr void set_max_K(ullong new_K) noexcept {
        K = new_K;
        return;
    }

    constexpr void set_vec (std::vector<T> new_vec) noexcept {
        vec = new_vec;
        return;
    }

    constexpr std::vector<T> return_vec() const noexcept{
        return vec;
    }

    void print_matrix() const noexcept {
        if (K == 1) {
            for (auto i = 0; i < I; ++i) {
                for (auto j = 0; j < J; ++j) {
                    std::cout<<vec[i + (I * j)]<<"\t";
                }
                std::cout<<"\n";
            }
            std::cout<<"\n";
            return;
        } else {
        std::cout<<"Error in print_matrix(): matrix must be 1D or 2D\n";
        return;
        }
    }


private:
    ullong I = 1;
    ullong J = 1;
    ullong K = 1;
    std::vector<T> vec;
};

template<typename T>
Flat_vec<T> matrix_mult(Flat_vec<T> const& a ,Flat_vec<T> const& b) noexcept {
    if (a.max_K() == 1 && b.max_K() == 1) {
        if (a.max_J() == b.max_I()) {
            Flat_vec<T> result {a.max_I(), b.max_J()};

            for (auto f = 0; f < a.max_I(); ++f) {
                for (auto g = 0; g < b.max_J(); ++g) {
                    for (auto i = 0; i < a.max_J(); ++i) {
                        result(f,g) += a(f, i) * b(i,g);
                    }
                }
            }

            return result;
        } else {
            std::cout<<"Error, incompatible matrix dimensions, returning empty matrix\n";
            return Flat_vec<T> {1,1,1};
        }
    } else {
        std::cout<<"Error: matrix_mult() undefined for 3 dimensional case, returning empty matrix\n";
        return Flat_vec<T> {1,1,1};
    }
}

#endif // FLAT_VEC_H_INCLUDED
