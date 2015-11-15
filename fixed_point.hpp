#pragma once

#include <utility>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>

#include "derivative.hpp"

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
    }

    void test_fp(double x0, double abstol)
    {
        const auto g1 = [](double x) -> double {
            return (x * x + 2) / 3;
        };

        // will result in NaNs because the value under the square root
        // will become negative at some point.
        const auto g2 = [](double x) -> double {
            return std::sqrt(3 * x - 2);
        };

        const auto g3 = [](double x) -> double {
            return 3 - (2 / x);
        };

        const auto g4 = [](double x) -> double {
            return (x * x  - 2) / (2 * x - 3);
        };

        std::vector<std::function<double(double)>> vec {g1, g2, g3, g4};
        std::vector<std::string> filenames {"fpg1.txt", "fpg2.txt", "fpg3.txt", "fpg4.txt"};
        std::vector<std::string> funcnames {"g1", "g2", "g3", "g4"};

        for (auto i = 0; i < vec.size(); ++i) {
            fp::test_fixed_point(vec[i], x0, abstol, funcnames[i], filenames[i]);
        }
    }

    template<typename Func, typename Float>
    std::tuple<std::vector<Float>, int, std::vector<Float>> newton_method(const Func& f, Float x0, Float abstol)
    {
        // vector to store intermediate x_i values
        std::vector<Float> xvec;
        std::vector<Float> rvec; // vector of rate approximations

        // do the first iteration outside the while loop
        auto x1 = x0 - (f(x0) / derivative(f, x0));
        xvec.push_back(x1);
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
            x_i = x_iminus1 - (f(x_iminus1) / derivative(f, x_iminus1));
            xvec.push_back(x_i);
            x_iplus1 = x_i - (f(x_i) / derivative(f, x_i));
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

    template<typename Func, typename Float>
    void test_newton_method(const Func& f, Float x0, Float abstol,
                            const std::string& funcname = "f", const std::string& filename = "test_f.txt")
    {
        const int nameWidth     = 24;
        const int numWidth      = 25;

        std::ofstream file;
        file.open(filename.c_str(), std::ios::out | std::ios::app);
        file << std::scientific << std::setprecision(15);
        file << "Getting the roots of '" << funcname
        << "' given x_0 = " << x0 << " and abstol = " << abstol << std::endl;

        auto tuple = newton_method(f, x0, abstol);

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
    }

    namespace newton {
        // we can't use functors with the derivative function
        // so we have to stick to traditional functions
        double f1(double x)
        {
            return x * x - 3 * x + 2; // q1
        }

        double f2(double x)
        {
            return x * x * x - 2 * x - 5; // q3a
        }

        double f3(double x)
        {
            return std::exp(-x) - x; // q3b
        }

        double f4(double x)
        {
            return std::sin(x) * x - 1; // q3c
        }

        double f5(double x)
        {
            return x * x * x - 3 * x * x + 3 * x - 1; // q3d
        }
    }

    void test_newton(double abstol)
    {
        using namespace newton; // for functions f1 to f5

        std::vector<double (*)(double)> vec {f1, f2, f3, f4, f5};
        std::vector<std::string> filenames {"newtonf1.txt", "newtonf2.txt",
                                            "newtonf3.txt", "newtonf4.txt", "newtonf5.txt"};
        std::vector<std::string> funcnames {"f1", "f2", "f3", "f4", "f5"};
        std::vector<double> xnaughts {2.1, 2.5, 0.6, 0.9, 0.5};

        for (auto i = 0; i < vec.size(); ++i) {
            test_newton_method(*vec[i], xnaughts[i], abstol, funcnames[i], filenames[i]);
        }
    }

    template<typename Func, typename Float>
    std::tuple<std::vector<Float>,
            int,
            std::vector<Float>>
    secant_method(const Func& f, Float x0, Float x1, Float abstol)
    {
        // vector to store intermediate x_i values
        std::vector<Float> xvec;
        std::vector<Float> rvec; // vector of rate approximations

        // do the first iteration outside the while loop
        xvec.push_back(x0);
        xvec.push_back(x1);
        auto x2 = x1 - f(x1) * ((x1 - x0) / (f(x1) - f(x0)));
        xvec.push_back(x2);
        Float currtol = static_cast<Float>(std::abs(xvec[2] - x1));

        // number of iterations : we have already performed one
        auto i = 3;
        Float rate = NAN;
        rvec.push_back(rate);

        Float x_iminus2;
        Float x_iminus1;
        Float x_i;
        Float x_iplus1;
        Float deltaxp1;
        while (currtol > abstol)
        {
            x_iminus2 = xvec[i - 2];
            x_iminus1 = xvec[i - 1];
            x_i = x_iminus1 - f(x_iminus1) * ((x_iminus1 - x_iminus2) / (f(x_iminus1) - f(x_iminus2)));
            xvec.push_back(x_i);
            x_iplus1 = x_i - f(x_i) * ((x_i - x_iminus1) / (f(x_i) - f(x_iminus1)));
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

    template<typename Func, typename Float>
    void test_secant_method(const Func& f, Float x0, Float x1, Float abstol,
                            const std::string& funcname = "f", const std::string& filename = "test_f.txt")
    {
        const int nameWidth     = 24;
        const int numWidth      = 25;

        std::ofstream file;
        file.open(filename.c_str(), std::ios::out | std::ios::app);
        file << std::scientific << std::setprecision(15);
        file << "Getting the roots of '" << funcname
        << "' given x_0 = " << x0 << ", x_1 = " << x1 << " and abstol = " << abstol << std::endl;

        auto tuple = secant_method(f, x0, x1, abstol);

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
    }

    void test_secant(double abstol)
    {
        using namespace newton; // using the newton functions

        std::vector<double (*)(double)> vec {f1, f2, f3, f4, f5};
        std::vector<std::string> filenames {"secantf1.txt", "secantf2.txt",
                                            "secantf3.txt", "secantf4.txt", "secantf5.txt"};
        std::vector<std::string> funcnames {"f1", "f2", "f3", "f4", "f5"};
        std::vector<double> xnaughts {2.5, -15, -1, 0, -4}; // xnaughts for secant
        std::vector<double> xones {2.1, -12.5, -0.5, 0.3, -3}; // xones for secant

        for (auto i = 0; i < vec.size(); ++i) {
            test_secant_method(*vec[i], xnaughts[i], xones[i], abstol, funcnames[i], filenames[i]);
        }
    }

    template<typename Float, typename = typename std::enable_if<std::is_arithmetic<Float>::value>::type>
    Float midpoint(Float a, Float b)
    {
        return a + (b - a) / 2;
    }

    template<typename Float, typename = typename std::enable_if<std::is_arithmetic<Float>::value>::type>
    int sign(Float a)
    {
        return a > 0 ? 1 : -1;
    }

    template<typename Func, typename Float>
    std::tuple<std::vector<Float>,
            int,
            std::vector<Float>>
    bisection_method(const Func& f, Float a, Float b, Float abstol, int numiters)
    {
        std::vector<Float> xvec; // to store intermediate root approximations
        std::vector<Float> rvec; // to store rate approximations
        auto n = 1;

        Float prevc;
        Float c;
        Float nextc;
        Float l = a;
        Float u = b;
        Float rate;
        Float currtol;
        while (n <= numiters)
        {
            prevc = c; // will be garbage in first iteration
            c = midpoint(l, u);
            xvec.push_back(u);
            currtol = std::abs(c - prevc);
            if (currtol < abstol) {
                break;
            }
            if (sign(f(c)) == sign(f(l))) {
                l = c;
            } else {
                u = c;
            }
            // get nextc to calculate rate
            nextc = midpoint(l, u);
            // will probably be garbage in the first iteration
            if (n >= 2) {
                rate = std::log(std::abs(nextc - c) / currtol) / std::log(currtol / std::abs(prevc - xvec[n - 2]));
                rvec.push_back(rate);
            }
            n++;
        }

        return std::make_tuple(xvec, n, rvec);
    }

    template<typename Func, typename Float>
    void test_bisection_method(const Func& f, Float a, Float b, Float abstol, int numiters,
                               const std::string& funcname = "f", const std::string& filename = "test_f.txt")
    {
        const int nameWidth     = 24;
        const int numWidth      = 25;

        std::ofstream file;
        file.open(filename.c_str(), std::ios::out | std::ios::app);
        file << std::scientific << std::setprecision(15);
        file << "Getting the roots of '" << funcname
        << "' given a = " << a << ", b = " << b << ", numiters = " << numiters
            << " and abstol = " << abstol << std::endl;

        auto tuple = bisection_method(f, a, b, abstol, numiters);

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
    }

    void test_bisection(double abstol, int numiters)
    {
        using namespace newton; // using the newton functions

        std::vector<double (*)(double)> vec {f1, f2, f3, f4, f5};
        std::vector<std::string> filenames {"bisectionf1.txt", "bisectionf2.txt",
                                            "bisectionf3.txt", "bisectionf4.txt", "bisectionf5.txt"};
        std::vector<std::string> funcnames {"f1", "f2", "f3", "f4", "f5"};
        // interval data
        std::vector<double> as {1.5, 1, -1, 0, 0.5}; // lower bound limit for bisection
        std::vector<double> bs {2.5, 3, 2, 2, 1.5}; // upper bound limit for bisection

        for (auto i = 0; i < vec.size(); ++i) {
            test_bisection_method(*vec[i], as[i], bs[i], abstol, numiters, funcnames[i], filenames[i]);
        }
    }
}