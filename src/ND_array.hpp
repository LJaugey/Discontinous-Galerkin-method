#ifndef NDARRAY_HPP
#define NDARRAY_HPP

#include <cstddef>
#include <iostream>
#include <type_traits>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Array_Expression.hpp"
#include "Unary_Expression.hpp"
#include "Binary_Expression.hpp"



namespace ND {

template <size_t firstDim, size_t... RestDims>
class Array : public Array_Expression<ND::Array<firstDim, RestDims...>>
{
    template<size_t dim1, size_t... dimN>
    friend class Array;
public:

    static constexpr size_t N = sizeof...(RestDims) + 1;
    static constexpr size_t length = firstDim * (RestDims * ...);
    static constexpr size_t Dims[N] = {firstDim, RestDims...};

    typedef typename base_traits<Array>::terminal_type terminal_type;
    typedef typename base_traits<Array>::terminal_sub_type terminal_sub_type;

protected:

    double* data_;
    bool is_original;


public:

    inline const double get_element(size_t i) const     {   return data_[i];    }

    // Base constructor
    Array()
    {
        data_ = new double[length];

        is_original =  true;
    }
    Array(double val)
    {
        data_ = new double[length];

        std::fill_n(data_,length, val);

        is_original =  true;
    }

    // copy constructor
    Array(const ND::Array<firstDim, RestDims...>& other)
    : is_original(true)
    {
        data_ = new double[length];

        std::copy(other.data_, other.data_ + length, data_);
    }

    // Constructor from pointer
    Array(double* p, bool is_or = false)
    {
        data_ = p;
        is_original = is_or;
    }
    
    // copy assigment operator
    const ND::Array<firstDim, RestDims...>& operator=(const ND::Array<firstDim, RestDims...>& other)
    {
        std::copy(other.data_, other.data_ + length, data_);
        return *this;
    }
    const ND::Array<firstDim, RestDims...>& operator=(double val)
    {
        std::fill_n(data_,length, val);
        return *this;
    }
    // constructor from N-1 dimensional array
    Array(const ND::Array<RestDims...>& slice)
    {
        data_ = new double[length];
        for(size_t i = 0; i<firstDim; i++)
        {
            std::copy(slice.data_, slice.data_ + slice.length, data_+ i*(RestDims*...));
        }
    }

    // construct from Array_expressions
    // Shift can be used if N<expr.N (e.g. operator[] on expressions)
    template <typename E>
    Array(const Array_Expression<E>& expr, size_t shift = 0)
    : is_original(true)
    {
        data_ = new double[length];
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; ++i)
            {
                data_[i] = expr.get_element(shift + i);
            }
        }
        else
        {
            for (size_t i = 0; i < length; ++i)
            {
                data_[i] = expr.get_element(shift + i);
            }
        }
    }
    template <typename E>
    const Array& operator=(const Array_Expression<E>& expr)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; ++i)
            {
                data_[i] = expr.get_element(i);
            }
        }
        else
        {
            for (size_t i = 0; i < length; ++i)
            {
                data_[i] = expr.get_element(i);
            }
        }

        return *this;
    }

    // destructor
    ~Array() {  if(is_original)    delete[] data_;   }



    // access element
    template <typename... ind_type>
    inline double& operator()(ind_type... indices)
    {
        size_t offset = 0;
        size_t temp[N] = {static_cast<size_t>(indices)...};

        for (size_t i = 0; i < std::min(N, sizeof...(ind_type)); i++) {
            offset = offset * Dims[i] + temp[i];
        }
        return data_[offset];
    }
    template <typename... ind_type>
    inline const double operator()(ind_type... indices) const
    {
        size_t offset = 0;
        size_t temp[N] = {static_cast<size_t>(indices)...};

        for (size_t i = 0; i < std::min(N, sizeof...(ind_type)); i++) {
            offset = offset * Dims[i] + temp[i];
        }
        return data_[offset];
    }


    // access element
    inline ND::Array<RestDims...> operator[](size_t index)
    {
        return ND::Array<RestDims...>(data_ + index * (RestDims * ...), false);  // Guaranteed copy elision
    }
    inline const ND::Array<RestDims...> operator[](size_t index) const
    {
        return ND::Array<RestDims...>(data_ + index * (RestDims * ...), false);  // Guaranteed copy elision
    }


    const ND::Array<firstDim, RestDims...>& fill(double val)
    {
        std::fill_n(data_,length, val);
        return *this;
    }
    const ND::Array<firstDim, RestDims...>& fill(const ND::Array<firstDim, RestDims...>& other)
    {
        if(data_ != other.data_)    std::copy(other.data_, other.data_ + length, data_);
        return *this;
    }
    template<class E>
    const ND::Array<firstDim, RestDims...>& fill(const Array_Expression<E>& expr)
    {
        for (size_t i = 0; i < length; ++i)
        {
            data_[i] = expr.get_element(i);
        }

        return *this;
    }


    inline const size_t size(const size_t index = 0) const
    {
        return Dims[index];
    }
    




    // Arithmetic operations

    // += operator
    template<class E>
    requires std::is_same_v<typename E::terminal_type, terminal_type>
    const ND::Array<firstDim, RestDims...>& operator+=(const Array_Expression<E>& expr)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] += expr.get_element(i);
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] += expr.get_element(i);
            }
        }
        return *this;
    }
    // scalar
    const ND::Array<firstDim, RestDims...>& operator+=(double scalar)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] += scalar;
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] += scalar;
            }
        }
        return *this;
    }

    // -= operator
    template<class E>
    requires std::is_same_v<typename E::terminal_type, terminal_type>
    const ND::Array<firstDim, RestDims...>& operator-=(const Array_Expression<E>& expr)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] -= expr.get_element(i);
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] -= expr.get_element(i);
            }
        }
        return *this;
    }
    // scalar
    const ND::Array<firstDim, RestDims...>& operator-=(double scalar)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] -= scalar;
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] -= scalar;
            }
        }
        return *this;
    }

    // *= operator
    template<class E>
    requires std::is_same_v<typename E::terminal_type, terminal_type>
    const ND::Array<firstDim, RestDims...>& operator*=(const Array_Expression<E>& expr)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] *= expr.get_element(i);
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] *= expr.get_element(i);
            }
        }
        return *this;
    }
    // scalar
    const ND::Array<firstDim, RestDims...>& operator*=(double scalar)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] *= scalar;
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] *= scalar;
            }
        }
        return *this;
    }

    // /= operator
    template<class E>
    requires std::is_same_v<typename E::terminal_type, terminal_type>
    const ND::Array<firstDim, RestDims...>& operator/=(const Array_Expression<E>& expr)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] /= expr.get_element(i);
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] /= expr.get_element(i);
            }
        }
        return *this;
    }
    // scalar
    const ND::Array<firstDim, RestDims...>& operator/=(double scalar)
    {
        double inv_scal = 1.0/scalar;

        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] *= inv_scal;
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] *= inv_scal;
            }
        }
        return *this;
    }
};
template<size_t firstDim, size_t... RestDims>
struct base_traits<ND::Array<firstDim, RestDims...>>
{
    typedef ND::Array<firstDim, RestDims...> terminal_type;
    typedef ND::Array<RestDims...> terminal_sub_type;
};




// ostream
template <size_t firstDim, size_t... RestDims>
std::ostream& operator<<(std::ostream& output, const ND::Array<firstDim, RestDims...>& other)
{
    if (sizeof...(RestDims) > 0)
    {
        for(size_t i = 0; i<firstDim; i++)
            output<<other[i]<<std::endl;
    }
    else
    {
        for(size_t i = 0; i<firstDim; i++)
            output<<other[i]<<"\t";
    }

    return output;
}



template <size_t Dim>
class Array<Dim> : public Array_Expression<ND::Array<Dim>>
{
    template<size_t dim1, size_t... dimN>
    friend class Array;
public:

    static constexpr size_t N = 1;
    static constexpr size_t length = Dim;
    static constexpr size_t Dims[N] = {Dim};

    typedef typename base_traits<Array>::terminal_type terminal_type;
    typedef typename base_traits<Array>::terminal_sub_type terminal_sub_type;

protected:

    double* data_;
    bool is_original;

    
public:

    inline const double get_element(size_t i) const     {   return data_[i];    }

    // Base constructor
    Array()
    {
        data_ = new double[length];

        is_original =  true;
    }
    Array(double val)
    {
        data_ = new double[length];

        std::fill_n(data_,length, val);

        is_original =  true;
    }

    // copy constructor
    Array(const ND::Array<Dim>& other)
    : is_original(true)
    {
        data_ = new double[length];

        std::copy(other.data_, other.data_ + length, data_);
    }
    
    // Constructor from pointer
    Array(double* p, bool is_or = false)
    {
        data_ = p;
        is_original = is_or;
    }
    
    // copy assigment operator
    const ND::Array<Dim>& operator=(const ND::Array<Dim>& other)
    {
        std::copy(other.data_, other.data_ + length, data_);
        return *this;
    }
    const ND::Array<Dim>& operator=(double val)
    {
        std::fill_n(data_,length, val);
        return *this;
    }
    
    // construct from Array_expressions
    // Shift can be used if N<expr.N (e.g. operator[] on expressions)
    template <typename E>
    Array(const Array_Expression<E>& expr, size_t shift = 0)
    : is_original(true)
    {
        data_ = new double[length];
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; ++i)
            {
                data_[i] = expr.get_element(shift + i);
            }
        }
        else
        {
            for (size_t i = 0; i < length; ++i)
            {
                data_[i] = expr.get_element(shift + i);
            }
        }
    }
    template <typename E>
    const Array& operator=(const Array_Expression<E>& expr)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; ++i)
            {
                data_[i] = expr.get_element(i);
            }
        }
        else
        {
            for (size_t i = 0; i < length; ++i)
            {
                data_[i] = expr.get_element(i);
            }
        }

        return *this;
    }

    // destructor
    ~Array() {  if(is_original)    delete[] data_;   }



    // access element
    inline double& operator()(size_t index)                { return data_[index]; }
    inline const double operator()(size_t index) const     { return data_[index]; }
    inline double& operator[](size_t index)                { return data_[index]; }
    inline const double operator[](size_t index) const     { return data_[index]; }





    const ND::Array<Dim>& fill(double val)
    {
        std::fill_n(data_,length, val);
        return *this;
    }
    const ND::Array<Dim>& fill(const ND::Array<Dim>& other)
    {
        if(data_ != other.data_)    std::copy(other.data_, other.data_ + length, data_);
        return *this;
    }
    template<class E>
    const ND::Array<Dim>& fill(const Array_Expression<E>& expr)
    {
        for (size_t i = 0; i < length; ++i)
        {
            data_[i] = expr.get_element(i);
        }

        return *this;
    }


    inline const size_t size(const size_t index = 0) const
    {
        return Dims[index];
    }






    // Arithmetic operations


    // += operator
    const ND::Array<Dim>& operator+=(const ND::Array<Dim>& other)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] += other.data_[i];
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] += other.data_[i];
            }
        }
        return *this;
    }
    // scalar
    const ND::Array<Dim>& operator+=(double scalar)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] += scalar;
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] += scalar;
            }
        }
        return *this;
    }

    // -= operator
    const ND::Array<Dim>& operator-=(const ND::Array<Dim>& other)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] -= other.data_[i];
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] -= other.data_[i];
            }
        }
        return *this;
    }
    // scalar
    const ND::Array<Dim>& operator-=(double scalar)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] -= scalar;
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] -= scalar;
            }
        }
        return *this;
    }

    // *= operator
    const ND::Array<Dim>& operator*=(const ND::Array<Dim>& other)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] *= other.data_[i];
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] *= other.data_[i];
            }
        }
        return *this;
    }
    // scalar
    const ND::Array<Dim>& operator*=(double scalar)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] *= scalar;
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] *= scalar;
            }
        }
        return *this;
    }

    // /= operator
    const ND::Array<Dim>& operator/=(const ND::Array<Dim>& other)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] /= other.data_[i];
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] /= other.data_[i];
            }
        }
        return *this;
    }
    // scalar
    const ND::Array<Dim>& operator/=(double scalar)
    {
        double inv_scal = 1.0/scalar;

        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; i++) {
                data_[i] *= inv_scal;
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++) {
                data_[i] *= inv_scal;
            }
        }
        return *this;
    }
};
template<size_t Dim>
struct base_traits<ND::Array<Dim>>
{
    typedef ND::Array<Dim> terminal_type;
    typedef double terminal_sub_type;
};

}

#endif
