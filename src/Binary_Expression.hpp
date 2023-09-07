#ifndef BINARY_EXPRESSION_HPP
#define BINARY_EXPRESSION_HPP

#include <cstddef>
#include <type_traits>

#include "Array_Expression.hpp"
#include "Unary_Expression.hpp"

namespace ND {

template <class E1, class OP, class E2>
class Binary_Op : public Array_Expression<Binary_Op<E1,OP,E2>>
{
    const E1& arg1;
    const E2& arg2;

public:
    
    typedef typename base_traits<Binary_Op>::terminal_type terminal_type;
    typedef typename base_traits<Binary_Op>::terminal_sub_type terminal_sub_type;
    
    Binary_Op(const E1& a_1, const E2& a_2)
    :arg1(a_1),arg2(a_2)
    {}

    inline const double get_element(size_t i) const
    {
        if constexpr(std::is_arithmetic_v<E1>)
        {
            return OP::apply(arg1,arg2.get_element(i));
        }
        else if constexpr(std::is_arithmetic_v<E2>)
        {
            return OP::apply(arg1.get_element(i),arg2);
        }
        else
        {
            return OP::apply(arg1.get_element(i),arg2.get_element(i));
        }
    }

    template <typename... ind_type>
    inline const double operator()(ind_type... indices) const
    {
        if constexpr(std::is_arithmetic_v<E1>)
        {
            return OP::apply(arg1,arg2(indices...));
        }
        else if constexpr(std::is_arithmetic_v<E2>)
        {
            return OP::apply(arg1(indices...),arg2);
        }
        else
        {
            return OP::apply(arg1(indices...),arg2(indices...));
        }
    }
};
template <class E1, class OP, class E2>
struct base_traits<Binary_Op<E1,OP,E2>>
{
    typedef typename std::conditional<  std::is_arithmetic_v<E1>,
                                        base_traits<E2>,
                                        base_traits<E1>
                                        >::type::terminal_type   terminal_type;

    typedef typename terminal_type::terminal_sub_type terminal_sub_type;
};


struct Array_add
{
    static inline const double apply(const double u, const double v)   {   return u + v;  }
};
template <class LHS, class RHS>
static inline const Binary_Op<LHS,Array_add,RHS> operator+(const LHS& lhs, const RHS& rhs)
{
    return Binary_Op<LHS,Array_add,RHS>(lhs,rhs);
}

struct Array_sub
{
    static inline const double apply(const double u, const double v)  {   return u - v;  }
};
template <class LHS, class RHS>
static inline const Binary_Op<LHS,Array_sub,RHS> operator-(const LHS& lhs, const RHS& rhs)
{
    return Binary_Op<LHS,Array_sub,RHS>(lhs,rhs);
}

struct Array_mul
{
    static inline const double apply(const double u, const double v)  {   return u * v;  }
};
template <class LHS, class RHS>
static inline const Binary_Op<LHS,Array_mul,RHS> operator*(const LHS& lhs, const RHS& rhs)
{
    return Binary_Op<LHS,Array_mul,RHS>(lhs,rhs);
}

struct Array_div
{
    static inline const double apply(const double u, const double v)  {   return u / v;  }
};
template <class LHS, class RHS>
static inline const Binary_Op<LHS,Array_div,RHS> operator/(const LHS& lhs, const RHS& rhs)
{
    return Binary_Op<LHS,Array_div,RHS>(lhs,rhs);
}

}

#endif