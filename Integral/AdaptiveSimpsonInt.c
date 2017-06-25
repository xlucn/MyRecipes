/** @file AdaptiveSimpsonInt.c */
#include <math.h>
#include "NR.h"

/**
 * @brief Adaptive Simpson intergral, using recursion.
 * @param f the integrand function
 * @param a the start point
 * @param b the end point
 * @param TOL the tolerance of precesion
 * @returns the integral
 */
double AdaptiveSimpsonInt(double (*f)(double), double a, double b, double TOL)
{
    double Stotal = SimpsonInt(f, a, b);
    double Sleft  = SimpsonInt(f, a, (a + b) / 2);
    double Sright = SimpsonInt(f, (a + b) / 2, b);

    if (fabs(Stotal - Sleft - Sright) > 10 * TOL)
    {
        Sleft = AdaptiveSimpsonInt(f, a, (a + b) / 2, TOL / 2);
        Sright = AdaptiveSimpsonInt(f, (a + b) / 2, b, TOL / 2);
    }

    return Sleft + Sright;
}
