/** @file TrapezoidalInt.c */
#include "NR.h"

/**
 * @brief Trapezoidal integration
 * @param f the integrand function
 * @param a the start point
 * @param b the end point
 * @returns the integral
 */
double TrapezoidalInt(double (*f)(double), double a, double b)
{
    return (b - a) / 2 * (f(a) + f(b));
}
