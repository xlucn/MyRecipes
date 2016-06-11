#include <stdlib.h>
#include <math.h>
#include "NR.h"
/**
 * @brief Trapezoidal integration
 */
double TrapezoidalInt(double(*f)(double), double a, double b)
{
    return (b - a) / 2 * (f(a) + f(b));
}
