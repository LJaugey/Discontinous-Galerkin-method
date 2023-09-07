#ifndef UNARY_EXPRESSION_HPP
#define UNARY_EXPRESSION_HPP

#include <cstddef>
#include <type_traits>
#include <algorithm>

#include "Array_Expression.hpp"

namespace ND {

template <class OP, class E>
class Unary_Op : public Array_Expression<Unary_Op<OP,E>>
{
    const E& arg;

public:

    typedef typename base_traits<Unary_Op>::terminal_type terminal_type;
    typedef typename base_traits<Unary_Op>::terminal_sub_type terminal_sub_type;

    Unary_Op(const E& a)
    :arg(a)
    {}

    inline const double get_element(size_t i) const
    {
        return OP::apply(arg.get_element(i));
    }

    template <typename... ind_type>
    inline const double operator()(ind_type... indices) const
    {
        return OP::apply(arg(indices...));
    }
};
template <class OP, class E>
struct base_traits<Unary_Op<OP,E>>
{
    typedef typename base_traits<E>::terminal_type terminal_type;

    typedef typename terminal_type::terminal_sub_type terminal_sub_type;
};


struct Array_opp
{
    static inline const double apply(const double u)    {   return -u;  }
};
template <class RHS>
static inline const Unary_Op<Array_opp,RHS> operator-(const RHS& rhs)
{
    return Unary_Op<Array_opp,RHS>(rhs);
}

struct Array_abs
{
    static inline const double apply(const double u)    {   return std::abs(u);  }
};
template <class RHS>
static inline const Unary_Op<Array_abs,RHS> abs(const RHS& rhs)
{
    return Unary_Op<Array_abs,RHS>(rhs);
}


// trigonometric functions

// sin, asin, sinh, asinh
struct Array_sin
{
    static inline const double apply(const double u)    {   return std::sin(u);  }
};
template <class RHS>
static inline const Unary_Op<Array_sin,RHS> sin(const RHS& rhs)
{
    return Unary_Op<Array_sin,RHS>(rhs);
}
struct Array_asin
{
    static inline const double apply(const double u)    {   return std::asin(u);  }
};
template <class RHS>
static inline const Unary_Op<Array_asin,RHS> asin(const RHS& rhs)
{
    return Unary_Op<Array_asin,RHS>(rhs);
}
struct Array_sinh
{
    static inline const double apply(const double u)    {   return std::sinh(u);  }
};
template <class RHS>
static inline const Unary_Op<Array_sinh,RHS> sinh(const RHS& rhs)
{
    return Unary_Op<Array_sinh,RHS>(rhs);
}
struct Array_asinh
{
    static inline const double apply(const double u)    {   return std::asinh(u);  }
};
template <class RHS>
static inline const Unary_Op<Array_asinh,RHS> asinh(const RHS& rhs)
{
    return Unary_Op<Array_asinh,RHS>(rhs);
}

// cos, acos, cosh, acosh
struct Array_cos
{
    static inline const double apply(const double u)    {   return std::cos(u);  }
};
template <class RHS>
static inline const Unary_Op<Array_cos,RHS> cos(const RHS& rhs)
{
    return Unary_Op<Array_cos,RHS>(rhs);
}
struct Array_acos
{
    static inline const double apply(const double u)    {   return std::acos(u);  }
};
template <class RHS>
static inline const Unary_Op<Array_acos,RHS> acos(const RHS& rhs)
{
    return Unary_Op<Array_acos,RHS>(rhs);
}
struct Array_cosh
{
    static inline const double apply(const double u)    {   return std::cosh(u);  }
};
template <class RHS>
static inline const Unary_Op<Array_cosh,RHS> cosh(const RHS& rhs)
{
    return Unary_Op<Array_cosh,RHS>(rhs);
}
struct Array_acosh
{
    static inline const double apply(const double u)    {   return std::acosh(u);  }
};
template <class RHS>
static inline const Unary_Op<Array_acosh,RHS> acosh(const RHS& rhs)
{
    return Unary_Op<Array_acosh,RHS>(rhs);
}


// tan, atan, tanh, atanh
struct Array_tan
{
    static inline const double apply(const double u)    {   return std::tan(u);  }
};
template <class RHS>
static inline const Unary_Op<Array_tan,RHS> tan(const RHS& rhs)
{
    return Unary_Op<Array_tan,RHS>(rhs);
}
struct Array_atan
{
    static inline const double apply(const double u)    {   return std::atan(u);  }
};
template <class RHS>
static inline const Unary_Op<Array_atan,RHS> atan(const RHS& rhs)
{
    return Unary_Op<Array_atan,RHS>(rhs);
}
struct Array_tanh
{
    static inline const double apply(const double u)    {   return std::tanh(u);  }
};
template <class RHS>
static inline const Unary_Op<Array_tanh,RHS> tanh(const RHS& rhs)
{
    return Unary_Op<Array_tanh,RHS>(rhs);
}
struct Array_atanh
{
    static inline const double apply(const double u)    {   return std::atanh(u);  }
};
template <class RHS>
static inline const Unary_Op<Array_atanh,RHS> atanh(const RHS& rhs)
{
    return Unary_Op<Array_atanh,RHS>(rhs);
}


// exp
struct Array_exp
{
    static inline const double apply(const double u)    {   return std::exp(u);  }
};
template <class RHS>
static inline const Unary_Op<Array_exp,RHS> exp(const RHS& rhs)
{
    return Unary_Op<Array_exp,RHS>(rhs);
}


// log
struct Array_log
{
    static inline const double apply(const double u)    {   return std::log(u);  }
};
template <class RHS>
static inline const Unary_Op<Array_log,RHS> log(const RHS& rhs)
{
    return Unary_Op<Array_log,RHS>(rhs);
}

struct Array_log10
{
    static inline const double apply(const double u)    {   return std::log10(u);  }
};
template <class RHS>
static inline const Unary_Op<Array_log10,RHS> log10(const RHS& rhs)
{
    return Unary_Op<Array_log10,RHS>(rhs);
}

struct Array_log2
{
    static inline const double apply(const double u)    {   return std::log2(u);  }
};
template <class RHS>
static inline const Unary_Op<Array_log2,RHS> log2(const RHS& rhs)
{
    return Unary_Op<Array_log2,RHS>(rhs);
}


// round, floor, ceil
struct Array_round
{
    static inline const double apply(const double u)    {   return std::round(u);  }
};
template <class RHS>
static inline const Unary_Op<Array_round,RHS> round(const RHS& rhs)
{
    return Unary_Op<Array_round,RHS>(rhs);
}

struct Array_floor
{
    static inline const double apply(const double u)    {   return std::floor(u);  }
};
template <class RHS>
static inline const Unary_Op<Array_floor,RHS> floor(const RHS& rhs)
{
    return Unary_Op<Array_floor,RHS>(rhs);
}

struct Array_ceil
{
    static inline const double apply(const double u)    {   return std::ceil(u);  }
};
template <class RHS>
static inline const Unary_Op<Array_ceil,RHS> ceil(const RHS& rhs)
{
    return Unary_Op<Array_ceil,RHS>(rhs);
}

}

#endif