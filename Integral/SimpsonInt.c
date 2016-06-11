#include <stdlib.h>
#include <math.h>
#include "NR.h"

/**
 * @brief Simpson integration
 */
double SimpsonInt(double(*f)(double), double a, double b)
{
    return (b - a) / 6 * (f(a) + 4 * f((a + b) / 2) + f(b));
}
