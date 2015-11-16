#include <iostream>

#include "fixed_point.hpp"

const auto abstol = 5e-10;
const auto numiters = 80;

int main() {
    fp::test_secant(abstol);
    //fp::test_newton(abstol);
    //fp::test_fp(2.1, abstol);
    //fp::test_bisection(abstol, numiters);
}