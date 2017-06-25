/** @file ClassicRungeKutta.c */
#include "NR.h"

/**
 * @brief Classic Runge-Kutta Method
 * @param f right function of ODE
 * @param a lower limit of interval
 * @param b upper limit of interval
 * @param y0 initial value of f
 * @param N number of subintervals
 * @return array of ys
 */
double *ClassicRungeKutta(double (*f)(double, double), double a, double b, double y0, int N)
{
    double k[4];
    double h = (b - a) / N;
    double x = a;
    double y = y0;
    double* result = newArray1d(N + 1);
    result[0] = y0;

    for(int i = 0; i < N; i++)
    {
        k[0] = f(x, y);
        k[1] = f(x + h / 2, y + h * k[0] / 2);
        k[2] = f(x + h / 2, y + h * k[1] / 2);
        k[3] = f(x + h, y + h * k[2]);
        y += h / 6 * (k[0] + 2 * k[1] + 2 * k[2] + k[3]);
        x += h;
        result[i + 1] = y;
    }

    return result;
}

