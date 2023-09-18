#ifndef HELPER_HPP
#define HELPER_HPP

#include <stddef.h>
#include <type_traits>


#define PAR_SIZE 1024


#ifdef _OPENMP
    #include <omp.h>
    #define STRINGIFY(a) #a
    #define OMP_FOR(n) _Pragma(STRINGIFY(omp parallel for simd if(n>PAR_SIZE)))
    #define OMP_FOR_add(n,var) _Pragma(STRINGIFY(omp parallel for simd reduction(+:var) if(n>PAR_SIZE)))
    #define OMP_FOR_min(n,var) _Pragma(STRINGIFY(omp parallel for simd reduction(min:var) if(n>PAR_SIZE)))
    #define OMP_FOR_max(n,var) _Pragma(STRINGIFY(omp parallel for simd reduction(max:var) if(n>PAR_SIZE)))
#else
    #define omp_get_thread_num() 0
    #define OMP_FOR(n)
    #define OMP_FOR_add(n,var)
    #define OMP_FOR_min(n,var)
    #define OMP_FOR_max(n,var)
#endif



namespace ND {


// Base traits
template <class E> 
struct base_traits
{
    typedef E value_type;
};



// Array expression
template <class E>
class Array_Expression;


template <class E>
struct is_Array_Expression
{
    static constexpr bool value = std::is_base_of<Array_Expression<E>,E>::value;
};



// Unary operation
template <class OP, class E>
requires (ND::is_Array_Expression<E>::value)
class Unary_Op;


template <class OP, class E>
struct base_traits<Unary_Op<OP,E>>
{
    typename base_traits<E>::value_type arg_value_type;
    typedef decltype(OP::apply(arg_value_type)) value_type;


    template<typename any_type>
    using generic_terminal_type = typename base_traits<E>::generic_terminal_type<any_type>;

    template<typename any_type>
    using generic_terminal_sub_type = typename base_traits<E>::generic_terminal_sub_type<any_type>;
    

    typedef generic_terminal_type<value_type> terminal_type;
    typedef generic_terminal_sub_type<value_type> terminal_sub_type;

};



// Binary operation
template <class E1, class OP, class E2>
requires (  (ND::is_Array_Expression<E1>::value || ND::is_Array_Expression<E1>::value) ||
            (std::is_convertible<typename base_traits<E1>::value_type, typename base_traits<E2>::value_type>::value ||
             std::is_convertible<typename base_traits<E2>::value_type, typename base_traits<E1>::value_type>::value))
class Binary_Op;


template <class E1, class OP, class E2>
struct base_traits<Binary_Op<E1,OP,E2>>
{
    typename base_traits<E1>::value_type arg1_value_type;
    typename base_traits<E2>::value_type arg2_value_type;

    typedef decltype(OP::apply(arg1_value_type,arg2_value_type)) value_type;


    template<typename any_type>
    using generic_terminal_type = typename std::conditional< ND::is_Array_Expression<E1>::value,
                                                    base_traits<E1>,
                                                    base_traits<E2>
                                                    >::type::generic_terminal_type<any_type>;

    template<typename any_type>
    using generic_terminal_sub_type = typename std::conditional< ND::is_Array_Expression<E1>::value,
                                                        base_traits<E1>,
                                                        base_traits<E2>
                                                        >::type::generic_terminal_sub_type<any_type>;


    typedef generic_terminal_type<value_type> terminal_type;
    typedef generic_terminal_sub_type<value_type> terminal_sub_type;
};



// ND Array
template<typename T, size_t firstDim, size_t... RestDims>
class Array;


template<typename T, size_t firstDim, size_t... RestDims>
struct base_traits<Array<T, firstDim, RestDims...>>
{
    typedef T value_type;


    template<typename any_type>
    using generic_terminal_type = Array<any_type, firstDim, RestDims...>;

    template<typename any_type>
    using generic_terminal_sub_type = Array<any_type, RestDims...>;


    typedef generic_terminal_type<value_type> terminal_type;
    typedef generic_terminal_sub_type<value_type> terminal_sub_type;
};


template<typename T, size_t Dim>
struct base_traits<Array<T, Dim>>
{
    typedef T value_type;


    template<typename any_type>
    using generic_terminal_type = Array<any_type, Dim>;

    template<typename any_type>
    using generic_terminal_sub_type = any_type;


    typedef generic_terminal_type<value_type> terminal_type;
    typedef generic_terminal_sub_type<value_type> terminal_sub_type;
};



// Mask Array
template <typename T, size_t firstDim, size_t... RestDims>
class Mask_Array;


}


#endif