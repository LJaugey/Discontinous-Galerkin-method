#ifndef ARRAY_EXPRESSION_HPP
#define ARRAY_EXPRESSION_HPP

#include <cstddef>
#include <iostream>
#include <type_traits>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

#define PAR_SIZE 1024



namespace ND {

template <class E> 
struct base_traits;

template<typename T>
requires std::is_arithmetic_v<T>
struct base_traits<T>
{
    typedef double terminal_type;
};


template <class E>
class Array_Expression
{
public:
    
    typedef typename base_traits<E>::terminal_type terminal_type;
    typedef typename base_traits<E>::terminal_sub_type terminal_sub_type;

    
    inline const double get_element(size_t i) const   {  return static_cast<const E&>(*this).get_element(i);   }

    template <typename... ind_type>
    inline const double operator()(ind_type... indices) const   {   return static_cast<const E&>(*this)(indices...);    }


    inline const terminal_type eval() const
    {
        return terminal_type(*this);    // Guaranteed copy elision
    }

    // operator[] only collapses sub-array
    inline const terminal_sub_type operator[](size_t index) const
    {
        if constexpr (terminal_type::N==1)
        {
            return get_element(index);
        }
        else
        {
            return terminal_sub_type(*this, index*terminal_sub_type::length);  // Guaranteed copy elision
        }
    }

    



    // min
    const double min() const
    {
        double res = get_element(0);
        
        if constexpr(terminal_type::length>PAR_SIZE)
        {
            #pragma omp parallel for reduction(min:res) if(omp_get_num_threads() == 1)
            for (size_t i = 1; i < terminal_type::length; i++)
                res = std::min(res,get_element(i));
        }
        else {
            for (size_t i = 1; i < terminal_type::length; i++)
                res = std::min(res,get_element(i));
        }
        return res;
    }
    // max
    const double max() const
    {
        double res = get_element(0);
        
        if constexpr(terminal_type::length>PAR_SIZE)
        {
            #pragma omp parallel for reduction(max:res) if(omp_get_num_threads() == 1)
            for (size_t i = 1; i < terminal_type::terminal_type::length; i++)
                res = std::max(res,get_element(i));
        }
        else {
            for (size_t i = 1; i < terminal_type::terminal_type::length; i++) {
                res = std::max(res,get_element(i));
            }
        }
        return res;
    }
    // sum
    const double sum() const
    {
        double res = get_element(0);
        
        if constexpr(terminal_type::terminal_type::length>PAR_SIZE)
        {
            #pragma omp parallel for reduction(+:res) if(omp_get_num_threads() == 1)
            for (size_t i = 1; i < terminal_type::terminal_type::length; i++)
                res += get_element(i);
        }
        else {
            for (size_t i = 1; i < terminal_type::terminal_type::length; i++) {
                res += get_element(i);
            }
        }
        return res;
    }
    // mean
    inline const double mean() const   {   return (this->sum())/terminal_type::length; }

    const double stdev() const
    {
        double res = 0.0;
        double m = this->mean();
        
        if constexpr(terminal_type::length>PAR_SIZE)
        {
            #pragma omp parallel for reduction(+:res) if(omp_get_num_threads() == 1)
            for (size_t i = 0; i < terminal_type::length; i++)
                res += get_element(i)*get_element(i);
        }
        else {
            for (size_t i = 0; i < terminal_type::length; i++)
                res += get_element(i)*get_element(i);
        }
        return sqrt(res/terminal_type::length - m*m);
    }

};


template<class E>
inline const double min(const Array_Expression<E>& expr)   {   return expr.min();  }
template<class E>
inline const double max(const Array_Expression<E>& expr)   {   return expr.max();  }
template<class E>
inline const double sum(const Array_Expression<E>& expr)   {   return expr.sum();  }
template<class E>
inline const double mean(const Array_Expression<E>& expr)  {   return expr.mean(); }
template<class E>
inline const double stdev(const Array_Expression<E>& expr) {   return expr.stdev();}

template<class E>
inline std::ostream& operator<<(std::ostream& output, const Array_Expression<E>& expr)
{
    return output<<expr.eval();
}


//comparison operators
template<class E1, class E2>
static const bool operator==(const Array_Expression<E1>& expr1, const Array_Expression<E2>& expr2)
{
    if constexpr(!std::is_same_v<typename E1::terminal_type, typename E2::terminal_type>)  return false;

    for(size_t i=0; i<E1::terminal_type::length; i++)
    {
        if(expr1.get_element(i) != expr2.get_element(i))  return false;
    }

    return true;
}

}

#endif