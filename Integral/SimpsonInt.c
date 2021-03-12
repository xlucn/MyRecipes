/**
 * @file SimpsonInt.c
 * @brief SimpsonInt.c
 */

#include "NR.h"

/**
 * @brief Sinpson integration
 * @param f the integrand function
 * @param a the start point
 * @param b the end point
 * @returns the integral
 */
double SimpsonInt(double (*f)(double), double a, double b)
{
    return (b - a) / 6 * (f(a) + 4 * f((a + b) / 2) + f(b));
}
