/**
 * @file SecentMethod.c
 * @brief SecentMethod.c
 */

#include <math.h>
#include "NR.h"
#include "constants.h"

/**
 * @brief Secent method to solve a equation f = 0
 * @param f the function to be solved
 * @param x0 the first point
 * @param x1 the second point
 * @param eps the allowed tolerance
 * @returns the root of f = 0
 */
double SecentMethod(double (*f)(double), double x0, double x1, double eps)
{
    double y0 = f(x0);
    double y1 = f(x1);
    double p;
    for(int i = 0; i < ITER_LIM; i++)
    {
        p = x1 - y1 * (x1 - x0) / (y1 - y0);
        if(fabs(p - x1) < eps)
        {
            return p;
        }
        x0 = x1;
        x1 = p;
        y0 = y1;
        y1= f(x1);
    }
    return NAN;
}
