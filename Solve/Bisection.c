/** @file Bisection.c */
#include <math.h>
#include "NR.h"

/**
 * @brief Bisection method to find the root of a equation
 * @param f a function
 * @param a lower limit of interval
 * @param b upper limit of interval
 * @param eps tolerance
 * @returns root of function between a and b
 */
double Bisection(double (*f)(double), double a, double b, double eps)
{
    double p;
    while(fabs(b - a) > eps)
    {
        p = (a + b) / 2;
        (f(a) * f(p) < 0) ? (b = p) : (a = p);
    }
    return p;
}
