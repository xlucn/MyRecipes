/**
 * @file TwoStageRungeKutta.c
 * Two stage Runge Kutta method
 */
#include <stdlib.h>
#include "NR.h"

/**
 * @brief Two stage Runge Kutta method
 * @param N step number
 * @param y0 initial values
 * @param f derivative function. dy/dt=f(t,y)
 * @param a lower limit of interval
 * @param b upper limit of interval
 * @param param parameter that determines value of c1 and c2
 * @returns
 */
static double* TwoStageRungeKutta(int N, double y0, double a, double b, double (*f)(double, double), double para)
{
    double c1 = 1 -  0.5 / para;
    double c2 = 1 - c1;
    double k1;
    double k2;
    double h = (b - a) / N;
    double x = a;
    double y = y0;
    double* result = malloc(N * sizeof(double));

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
 * @param N step number
 * @param y0 initial values
 * @param f derivative function. dy/dt=f(t,y)
 * @param a lower limit of interval
 * @param b upper limit of interval
 * @returns
 */
double* ImprovedEuler(double (*f)(double, double), double a, double b, double y0, int N)
{
    return TwoStageRungeKutta(N, y0, a, b, f, 1);
}

/**
 * @brief Midpoint method.
 * @param N step number
 * @param y0 initial values
 * @param f derivative function. dy/dt=f(t,y)
 * @param a lower limit of interval
 * @param b upper limit of interval
 * @returns
 */
double* MID(double (*f)(double, double), double a, double b, double y0, int N)
{
    return TwoStageRungeKutta(N, y0, a, b, f, 1 / 2);
}

/**
 * @brief Heun Method
 * @param N step number
 * @param y0 initial values
 * @param f derivative function. dy/dt=f(t,y)
 * @param a lower limit of interval
 * @param b upper limit of interval
 * @returns
 */
double* Heun(double (*f)(double, double), double a, double b, double y0, int N)
{
    return TwoStageRungeKutta(N, y0, a, b, f, 2 / 3);
}
