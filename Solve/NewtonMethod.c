#include <math.h>
#include "NR.h"

/**
 * @brief Newton method (or Newton-Raphson Method) to find the root of a equation.
 * df is the derivative function of f(x)
 */
double NewtonMethod(double(*f)(double), double df(double), double x, double eps)
{
    double p = x - f(x) / df(x);
    while(fabs(p - x) > eps)
    {
        x = p;
        p = x - f(x) / df(x);
    }
    return p;
}
