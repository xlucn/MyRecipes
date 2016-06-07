#include <math.h>
#include "NR.h"

/**
 * @brief Accelerate iteration of PicardIteration
 */
double SteffensenIteration(double(*g)(double), double x, double eps)
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
