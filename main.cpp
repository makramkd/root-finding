#include <iostream>

#include "fixed_point.hpp"


int main() {
    auto g1 = [](double x) -> double {
        return (x * x + 2) / 3;
    };

    auto g2 = [](double x) -> double {
        return std::sqrt(3 * x - 2);
    };

    auto g3 = [](double x) -> double {
        return 3 - (2 / x);
    };

    auto g4 = [](double x) -> double {
        return (x * x  - 2) / (2 * x - 3);
    };

    std::vector<std::function<double(double)>> vec {g1, g2, g3, g4};
    std::vector<std::string> filenames {"g1.txt", "g2.txt", "g3.txt", "g4.txt"};
    std::vector<std::string> funcnames {"g1", "g2", "g3", "g4"};
    const double x0 = 0;
    const double abstol = 5e-7;

    for (auto i = 0; i < vec.size(); ++i) {
        fp::test_fixed_point(vec[i], x0, abstol, funcnames[i], filenames[i]);
    }
}