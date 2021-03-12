/**
 * @file PicardIteration.c
 * @brief PicardIteration.c
 */

#include <math.h>
#include "NR.h"

/**
 * @brief Use integration to solve equations in the form of x = g(x)
 * @param g function in fixed point equation
 * @param x initial point
 * @param eps tolerance
 * @returns the root of the equation
 */
double PicardIteration(double (*g)(double), double x, double eps)
{
    double p = g(x);
    while(fabs(p - x) > eps)
    {
        x = p;
        p = g(x);
    }
    return p;
}

