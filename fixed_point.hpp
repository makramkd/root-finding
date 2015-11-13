#pragma once

#include <utility>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>

namespace fp {

    template<typename Func, typename Float>
    std::tuple<std::vector<Float>, int, std::vector<Float>> fixed_point(const Func& g, Float x0, Float abstol)
    {
        // vector to store intermediate x_i values
        std::vector<Float> xvec;
        std::vector<Float> rvec; // vector of rate approximations

        // do the first iteration outside the while loop
        xvec.push_back(g(x0));
        Float currtol = static_cast<Float>(std::abs(xvec[0] - x0));

        // number of iterations : we have already performed one
        auto i = 1;
        Float rate = NAN;
        rvec.push_back(rate);

        Float x_iminus1;
        Float x_i;
        Float x_iplus1;
        Float deltaxp1;
        while (currtol > abstol)
        {
            x_iminus1 = xvec[i - 1];
            x_i = g(x_iminus1);
            xvec.push_back(x_i); // xi = g(x_{i - 1})
            x_iplus1 = g(x_i);
            currtol = std::abs(x_i - x_iminus1);
            deltaxp1 = std::abs(x_iplus1 - x_i);

            if ((i - 2) > 0) {
                rate = std::log(deltaxp1 / currtol) / std::log(currtol / std::abs(x_iminus1 - xvec[i - 2]));
                rvec.push_back(rate);
            } else {
                rvec.push_back(NAN);
            }

            ++i;
        }

        return std::make_tuple(xvec, i, rvec);
    }

    template<typename T> void printElement(T t, const int& width, std::ostream& stream)
    {
        static const char separator = ' ';
        stream << std::left << std::setw(width) << std::setfill(separator) << t;
    }

    template<typename Func, typename Float>
    void test_fixed_point(const Func& g, Float x0, Float abstol,
                          const std::string& funcname = "g", const std::string& filename = "test_g.txt")
    {
        const int nameWidth     = 24;
        const int numWidth      = 25;

        std::ofstream file;
        file.open(filename.c_str(), std::ios::out | std::ios::app);
        file << std::scientific << std::setprecision(15);
        file << "Getting the fixed points of '" << funcname
            << "' given x_0 = " << x0 << " and abstol = " << abstol << std::endl;

        auto tuple = fixed_point(g, x0, abstol);

        auto xvec = std::get<0>(tuple);
        auto iters = std::get<1>(tuple);
        auto rvec = std::get<2>(tuple);

        printElement("i", nameWidth, file);
        printElement("x_i", nameWidth, file);
        printElement("|x_i - x_{i - 1}|", nameWidth, file);
        printElement("rate", nameWidth, file);
        file << '\n';

        for (auto i = 0; i < xvec.size(); ++i)
        {
            printElement(i, numWidth, file);
            printElement(xvec[i], numWidth, file);
            if ((i - 1) >= 0) {
                printElement(std::abs(xvec[i] - xvec[i - 1]), numWidth, file);
            }
            printElement(rvec[i], numWidth, file);
            file << '\n';
        }

        file << "END" << std::endl;
    };
}