#include <math.h>
#include "NR.h"

/**
 * @brief Muller method
 */
double MullerMethod(double(*f)(double), double x0, double x1, double x2, double eps)
{
    double b;
    double d;
    double e;
    double h;

    double f0 = f(x0);
    double f1 = f(x1);
    double f2 = f(x2);
    double f01 = (f1 - f0) / (x1 - x0);
    double f12 = (f2 - f1) / (x2 - x1);
    double f012  = (f12 - f01) / (x2 - x0);

    do{
        b = f12 + f012 * (x2 - x1);
        d = sqrt(b * b - 4 * f2 * f012);
        e = fabs(b - d) > fabs(b + d) ? b + d : b - d;
        h = -2 * f2 / e;
        x0 = x1;
        x1 = x2;
        x2 += h;
        f0 = f1;
        f1 = f2;
        f2 = f(x2);
        f01 = (f1 - f0) / (x1 - x0);
        f12 = (f2 - f1) / (x2 - x1);
        f012  = (f12 - f01) / (x2 - x0);
    }while(fabs(h) < eps);

    return x2;
}
