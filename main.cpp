#include <iostream>

#include "derivative.hpp"

int f(int, int, int)
{
    return 0;
}

double g(double*, char*, std::ostream&)
{
    return 0.0;
}

double f1(double x)
{
    return x * x;
}

double f2(double x)
{
    return 2 * x * x + 3 * x;
}

int main() {

    auto df1 = fp::derivative(f1, 1.);
    auto df2 = fp::derivative(f2, 1.0 / 2.0);

    std::cout << df1 << std::endl;
    std::cout << df2 << std::endl;
}