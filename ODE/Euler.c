/**
 * @file Euler.c
 * @brief Euler.c
 */

#include <stdlib.h>
#include "NR.h"

/**
 * @brief Euler method to solve initial value problem(IVP) of ODE
 * @param f derivative function. dy/dt=f(t,y)
 * @param a lower limit of interval
 * @param b upper limit of interval
 * @param y0 initial value
 * @param N number of subintervals
 * @return
 * @note you have to free the returning pointer after usage
 */
double* Euler(double (*f)(double, double), double a, double b, double y0, int N)
{
    double h = (b - a) / N;
    double x = a;
    double* y = malloc((N + 1) * sizeof(double));
    y[0] = y0;

    for(int i = 0; i < N; i++)
    {
        y[i + 1] += h * f(x, y[i]);
        x += h;
    }
    return y;
}
