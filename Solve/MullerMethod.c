/** @file MullerMethod.c */
#include <math.h>
#include "NR.h"
#include "constants.h"

/**
 * @brief solve a equation with Muller method
 * @param f the function in equation f = 0
 * @param x0 the first initial point
 * @param x1 the second initial point
 * @param x2 the third initial point
 * @param eps the tolerance allowed
 * @returns the root of the function, macro NAN will be returned if the
 *  limit of iteration #ITER_LIM is reached.
 */
double MullerMethod(double (*f)(double), double x0, double x1, double x2, double eps)
{
    double b, d, e, h;
    double f0, f1, f2, f01, f12, f012;
    f0 = f(x0);f1 = f(x1);

    for(int i = 0; i < ITER_LIM; i++)
    {
        f2 = f(x2);
        f01 = (f1 - f0) / (x1 - x0);
        f12 = (f2 - f1) / (x2 - x1);
        f012  = (f12 - f01) / (x2 - x0);
        b = f12 + f012 * (x2 - x1);
        d = sqrt(b * b - 4 * f2 * f012);
        e = fabs(b - d) < fabs(b + d) ? b + d : b - d;
        h = -2 * f2 / e;
        if(fabs(h) < eps)
        {
            return x2;
        }
        x0 = x1;
        f0 = f1;
        x1 = x2;
        f1 = f2;
        x2 += h;
    }

    return NAN;
}
