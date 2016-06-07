#include <math.h>
#include "NR.h"

/**
 * @brief Bisection method to find the root of a equation
 */
double Bisection(double(*f)(double), double a, double b, double eps)
{
    double p;
    while(fabs(b - a) > eps)
    {
        p = (a + b) / 2;
        (f(a) * f(p) < 0) ? (b = p) : (a = p);
    }
    return p;
}
