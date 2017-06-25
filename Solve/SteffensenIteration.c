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
    double y = g(x);
    double z = g(y);
    double p = x - (y - x) * (y - x) / (z - 2 * y + x);
    while(fabs(p - x) > eps)
    {
        x = p;
        y = g(x);
        z = g(y);
        p = x - (y - x) * (y - x) / (z - 2 * y + x);
    }
    return p;
}
