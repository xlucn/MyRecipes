#include <math.h>
#include "NR.h"

/**
 * @brief Use integration to solve equations in the form of x = g(x)
 */
double PicardIteration(double(*g)(double), double x, double eps)
{
    double p = g(x);
    while(fabs(p - x) > eps)
    {
        x = p;
        p = g(x);
    }
    return p;
}

