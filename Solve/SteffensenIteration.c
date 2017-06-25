/**
 * @file SteffensenIteration.c
 * @brief Accelerate iteration using PicardIteration
 */
#include <math.h>
#include "NR.h"
#include "constants.h"

/**
 * @brief Accelerate iteration using PicardIteration
 * @param g function as in x = g(x)
 * @param x the initial point
 * @param eps the allows tolerance
 * @returns the root of the function
 */
double SteffensenIteration(double (*g)(double), double x, double eps)
{
    double y, z, p;
    for(int i = 0; i < ITER_LIM; i++)
    {
        y = g(x);
        z = g(y);
        p = x - (y - x) * (y - x) / (z - 2 * y + x);
        if(fabs(p - x) < eps)
        {
            return p;
        }
        x = p;
    }
    return FLOAT_NAN;
}
