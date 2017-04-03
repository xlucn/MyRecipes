#include "NR.h"

/**
 * @brief Composite Simpson integration，N is the number of the subintervals
 */
double CompositeSimpsonInt(double(*f)(double), double a, double b, int N)
{
    double Integral = 0;
    double h = (b - a) / N / 2;

    for (int i = 0; i < N; i++)
    {
        Integral += f(a + 2 * i * h) * h / 3;
        Integral += f(a + (2 * i + 1) * h) * 4 * h / 3;
        Integral += f(a + (2 * i + 2) * h) * h / 3;
    }

    return Integral;
}
