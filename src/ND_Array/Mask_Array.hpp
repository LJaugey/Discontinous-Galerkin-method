#ifndef MASK_ARRAY_HPP
#define MASK_ARRAY_HPP

#include "helper.hpp"


namespace ND {


// Mask Array
template <typename T, size_t firstDim, size_t... RestDims>
class Mask_Array
{
    template<typename T_, size_t f_Dim, size_t... R_dims>
    friend class Array;

protected:

    static constexpr size_t length = firstDim * (RestDims * ...);
    using value_type = T;

    Array<T, firstDim, RestDims...>& arr_;
    const Array<bool, firstDim, RestDims...>& mask_;


    Mask_Array(Array<T, firstDim, RestDims...>& arr, const Array<bool, firstDim, RestDims...>& mask)
    : arr_(arr), mask_(mask)
    {}

public:
    template<class E>
    void operator=(const Array_Expression<E>& expr)
    {
        OMP_FOR(length)
        for (size_t i = 0; i < length; ++i)
            if(mask_.get_element(i))    arr_.data_[i] = expr.get_element(i);
    }
    void operator=(const value_type val)
    {
        OMP_FOR(length)
        for (size_t i = 0; i < length; ++i)
            if(mask_.get_element(i))    arr_.data_[i] = val;
    }
    template<class E>
    void operator+=(const Array_Expression<E>& expr)
    {
        OMP_FOR(length)
        for (size_t i = 0; i < length; ++i)
            if(mask_.get_element(i))    arr_.data_[i] += expr.get_element(i);
    }
    template<class E>
    void operator-=(const Array_Expression<E>& expr)
    {
        OMP_FOR(length)
        for (size_t i = 0; i < length; ++i)
            if(mask_.get_element(i))    arr_.data_[i] -= expr.get_element(i);
    }
    template<class E>
    void operator*=(const Array_Expression<E>& expr)
    {
        OMP_FOR(length)
        for (size_t i = 0; i < length; ++i)
            if(mask_.get_element(i))    arr_.data_[i] *= expr.get_element(i);
    }
    template<class E>
    void operator/=(const Array_Expression<E>& expr)
    {
        OMP_FOR(length)
        for (size_t i = 0; i < length; ++i)
            if(mask_.get_element(i))    arr_.data_[i] /= expr.get_element(i);
    }
};


template <typename T, size_t Dim>
class Mask_Array<T, Dim>
{
    template<typename T_, size_t f_Dim, size_t... R_dims>
    friend class Array;

protected:

    static constexpr size_t length = Dim;
    using value_type = T;

    Array<T, Dim>& arr_;
    const Array<bool, Dim>& mask_;


    Mask_Array(Array<T,Dim>& arr, const Array<bool,Dim>& mask)
    : arr_(arr), mask_(mask)
    {}

public:
    template<class E>
    void operator=(const Array_Expression<E>& expr)
    {
        OMP_FOR(length)
        for (size_t i = 0; i < length; ++i)
            if(mask_.get_element(i))    arr_.data_[i] = expr.get_element(i);
    }
    void operator=(const value_type val)
    {
        OMP_FOR(length)
        for (size_t i = 0; i < length; ++i)
            if(mask_.get_element(i))    arr_.data_[i] = val;
    }
    template<class E>
    void operator+=(const Array_Expression<E>& expr)
    {
        OMP_FOR(length)
        for (size_t i = 0; i < length; ++i)
            if(mask_.get_element(i))    arr_.data_[i] += expr.get_element(i);
    }
    template<class E>
    void operator-=(const Array_Expression<E>& expr)
    {
        OMP_FOR(length)
        for (size_t i = 0; i < length; ++i)
            if(mask_.get_element(i))    arr_.data_[i] -= expr.get_element(i);
    }
    template<class E>
    void operator*=(const Array_Expression<E>& expr)
    {
        OMP_FOR(length)
        for (size_t i = 0; i < length; ++i)
            if(mask_.get_element(i))    arr_.data_[i] *= expr.get_element(i);
    }
    template<class E>
    void operator/=(const Array_Expression<E>& expr)
    {
        OMP_FOR(length)
        for (size_t i = 0; i < length; ++i)
            if(mask_.get_element(i))    arr_.data_[i] /= expr.get_element(i);
    }
    
    ~Mask_Array()
    {}
};


}


#endif