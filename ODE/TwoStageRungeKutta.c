#include "NR.h"
#include "NRprivate.h"

/**
 * Two stage Runge Kutta method
 */
static double* TwoStageRungeKutta(int N, double y0, double a, double b, double(*f)(double, double), double para)
{
    double c1 = 1 -  0.5 / para;
    double c2 = 1 - c1;
    double k1;
    double k2;
    double h = (b - a) / N;
    double x = a;
    double y = y0;
    double* result = newArray1d(N);

    for(int i = 0; i < N; i++)
    {
        k1 = f(x, y);
        k2 = f(x + para * h, y + para * h * k1);
        y += h * (c1 * k1 + c2 * k2);
        x += h;
        result[i] = y;
    }
    return result;
}
/**
 * @brief Improved Euler Method
 */
double* ImprovedEuler(double(*f)(double, double), double a, double b, double y0, int N)
{
    return TwoStageRungeKutta(N, y0, a, b, f, 1);
}
/**
 * @brief Midpoint method.
 */
double* MID(double(*f)(double, double), double a, double b, double y0, int N)
{
    return TwoStageRungeKutta(N, y0, a, b, f, 1 / 2);
}
/**
 * @brief Heun Method
 */
double* Heun(double(*f)(double, double), double a, double b, double y0, int N)
{
    return TwoStageRungeKutta(N, y0, a, b, f, 2 / 3);
}
