/**
 * @file NewtonMethod.c
 * @brief NewtonMethod.c
 */

#include <math.h>
#include "NR.h"

/**
 * @brief Newton method (or Newton-Raphson Method) to find the root of a equation.
 * @param f a function
 * @param df derivative of f
 * @param x initial point
 * @param eps tolerance
 * @returns root of equation f = 0
 */
double NewtonMethod(double (*f)(double), double (*df)(double), double x, double eps)
{
    double p = x - f(x) / df(x);
    while(fabs(p - x) > eps)
    {
        x = p;
        p = x - f(x) / df(x);
    }
    return p;
}
