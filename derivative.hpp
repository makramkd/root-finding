#ifndef FP_DERIVATIVE_HPP
#define FP_DERIVATIVE_HPP

#include <type_traits>
#include <cmath>
#include <limits>

namespace fp {

    // we'd like a template to determine the arity
    // of a function/functor f. We can do this with
    // a metafunction.
    template<typename Func>
    struct arity;

    template<typename Ret, typename... Args>
    struct arity<Ret(Args...)>
    {
        static constexpr auto value = sizeof...(Args);
    };

    // we'd like a template to determine that the
    // return type of a function is an arithmetic type
    // this way we can use it in an enable_if
    template<typename Func>
    struct is_arithmetic_function;

    // use std::is_arithmetic so that we don't need
    // to do a lot of specializations
    template<typename Ret, typename... Args>
    struct is_arithmetic_function<Ret(Args...)>
    {
        static constexpr auto value = std::is_arithmetic<Ret>::value;
    };

    // dummy template
    template<typename Ret, typename... Args>
    struct arithmetic_args;

    // a template to determine whether the arguments
    // of a function are all arithmetic types
    template<typename Ret, typename Arg0, typename... Args>
    struct arithmetic_args<Ret(Arg0, Args...)>
    {
        static constexpr auto value = std::is_arithmetic<Arg0>::value && arithmetic_args<Ret(Args...)>::value;
    };

    template<typename Ret, typename Arg>
    struct arithmetic_args<Ret(Arg)>
    {
        static constexpr auto value = std::is_arithmetic<Arg>::value;
    };

    // get the return type of a function
    template<typename Func>
    struct return_type;

    template<typename Ret, typename... Args>
    struct return_type<Ret(Args...)>
    {
        using type = Ret;
    };

    // derivative routine:
    // check if the function has arity 1 and returns an arithmetic type
    // otherwise substitution failure will kick in
    template<typename F>
    typename std::enable_if<((arity<F>::value == 1) &&
            (is_arithmetic_function<F>::value)),
            typename return_type<F>::type>::type derivative(const F& func, typename return_type<F>::type x)
    {
        const auto h = std::sqrt(std::numeric_limits<typename return_type<F>::type>::epsilon());
        return (func(x + h) - func(x - h)) / (2 * h);
    };
}

#endif