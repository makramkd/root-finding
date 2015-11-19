#include <iostream>

#include "fixed_point.hpp"

const auto abstol = 5e-10;
const auto numiters = 80;

int main() {
    fp::test_newton_2(abstol);
}