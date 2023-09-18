#ifndef DYN_ARRAY_HPP
#define DYN_ARRAY_HPP

#include <iostream>
#include <type_traits>

#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Array_Expression.hpp"
#include "Unary_Expression.hpp"
#include "Binary_Expression.hpp"

#include <type_traits>
#include <typeinfo>
#ifndef _MSC_VER
#   include <cxxabi.h>
#endif
#include <memory>
#include <string>
#include <cstdlib>

template <class T>
std::string
type_name()
{
    typedef typename std::remove_reference<T>::type TR;
    std::unique_ptr<char, void(*)(void*)> own
           (
#ifndef _MSC_VER
                abi::__cxa_demangle(typeid(TR).name(), nullptr,
                                           nullptr, nullptr),
#else
                nullptr,
#endif
                std::free
           );
    std::string r = own != nullptr ? own.get() : typeid(TR).name();
    if (std::is_const<TR>::value)
        r += " const";
    if (std::is_volatile<TR>::value)
        r += " volatile";
    if (std::is_lvalue_reference<T>::value)
        r += "&";
    else if (std::is_rvalue_reference<T>::value)
        r += "&&";
    return r;
}


//namespace ND {
/*
template <typename T>
class Array : public Array_Expression<Array<T>>
{
    template<typename T_, size_t... Any_dims>
    friend class Array;

    
protected:
    
    size_t N;
    size_t length;
    size_t capacity;

    size_t* Dims;

    value_type* data_;
    bool is_original;

public:

    typedef typename base_traits<Array>::terminal_type terminal_type;
    typedef typename base_traits<Array>::terminal_sub_type terminal_sub_type;
    typedef typename base_traits<Array>::value_type  value_type;


    inline const value_type get_element(size_t i) const     {   return data_[i];    }


    template <typename... ind_type>
    Array(ind_type... indices)
    : is_original(true)
    {
        static_assert<size_t>(indices)...

        N = sizeof...(indices);
        Dims = {static_cast<size_t>(indices)...};
        length = (static_cast<size_t>(indices) * ...);
        capacity = length;

        data_ = new value_type[capacity];
    }

    // copy constructor
    Array(const Array<T>& other)
    : is_original(true)
    {
        N = other.N;
        std::copy(other.Dims, other.Dims + N, Dims);
        length = other.length;
        capacity = length;

        data_ = new value_type[other.capacity];

        std::copy(other.data_, other.data_ + length, data_);
    }

    // copy assigment operator
    const Array<T, firstDim, RestDims...>& operator=(const Array<T, firstDim, RestDims...>& other)
    {
        N = other.N;
        std::copy(other.Dims, other.Dims + N, Dims);
        length = other.length;
        capacity = length;

        data_ = new value_type[other.capacity];

        std::copy(other.data_, other.data_ + length, data_);

        return *this;
    }
    const Array<T, firstDim, RestDims...>& operator=(value_type val)
    {
        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < length; ++i)
            {
                data_[i] = val;
            }
        }
        else
        {
            std::fill_n(data_,length, val);
        }
        
        return *this;
    }
    // constructor from N-1 dimensional array
    Array(const Array<T, RestDims...>& slice)
    : is_original(true)
    {
        data_ = new value_type[length];

        if constexpr(length>PAR_SIZE)
        {
            #pragma omp parallel for if(omp_get_num_threads() == 1)
            for(size_t i = 0; i<length; i++)
            {
                data_[i] = slice.data_[i%(RestDims*...)];
            }
        }
        else
        {
            for(size_t i = 0; i<length; i++)
            {
                data_[i] = slice.data_[i%(RestDims*...)];
            }
        }
        
    }

    // construct from Array_expressions
    // Shift can be used if N<expr.N (e.g. operator[] on expressions)
    template <typename E>
    Array(const Array_Expression<E>& expr, size_t shift = 0)
    : is_original(true)
    {
        data_ = new value_type[length];
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

};
template<typename T>
struct base_traits<Array<T>>
{
    typedef Array<T> terminal_type;
    typedef Array<T> terminal_sub_type;
    typedef T value_type;
};
*/



/*
template <typename T, size_t N>
class Dynamic_Array : public Array_Expression<Dynamic_Array<T,N>>
{
public:
    Dynamic_Array()
    {
        std::cout<<"Type: "<<typeid(T).name()<<std::endl;
        std::cout<<"Rank: "<<N<<std::endl;
    }
};

template<typename T, size_t Rank>
struct base_traits<Dynamic_Array<T, Rank>>
{
    typedef T value_type;


    template<typename any_type>
    using generic_terminal_type = Dynamic_Array<any_type, Rank>;

    template<typename any_type>
    using generic_terminal_sub_type = Dynamic_Array<any_type, Rank-1>;


    typedef generic_terminal_type<value_type> terminal_type;
    typedef generic_terminal_sub_type<value_type> terminal_sub_type;
};


template<typename T>
using Array<T> = Dynamic_Array<T,0>;


template<typename T, size_t...D>
requires ((D == 0) && ...)
class jaj
{
    static constexpr size_t N = sizeof...(D);
public:
    jaj()
    {
        std::cout<<"Type: "<<typeid(T).name()<<std::endl;
        std::cout<<"Rank: "<<N<<std::endl;
    }
};
*/



template<typename T>
class jaj
{
public:
    T* data;
    int length;
    std::vector<int> dimensions;


    jaj(std::vector<int> dims) {
        length = 1;
        for (int dim : dims) {
            dimensions.push_back(dim);
            length *= dim;
        }
        data = new T[length];
    }
    jaj(std::vector<int> dims, T* data_P) {
        length = 1;
        for (int dim : dims) {
            dimensions.push_back(dim);
            length *= dim;
        }
        data = data_P;
    }

    // Access element at given indices
    T& operator()(std::vector<int> indices) {
        int index = 0;
        int multiplier = 1;
        for (int i = dimensions.size() - 1; i >= 0; i--) {
            index += indices[i] * multiplier;
            multiplier *= dimensions[i];
        }
        return data[index];
    }

    // Get dimensions of the array
    std::vector<int> getDimensions() const {
        return dimensions;
    }

    // Operator[] to return a slice or value
    auto operator[](int index) {
        if (dimensions.size() == 1) {
            return data[index];
        } else {
            std::vector<int> newDims(dimensions.begin() + 1, dimensions.end());
            int startIndex = index * (length / dimensions[0]);
            return jaj<T>(newDims, data + startIndex);
        }
    }
};

//}

#endif