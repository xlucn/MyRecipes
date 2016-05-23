#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "NR.h"
#include "Test.h"

/**
 * testing solving a equations
 */

static double f9(double x)
{
    return x * x * x - x - 1;
}

/**
 * @brief 
 * @returns 
 * 
 * 
 */
int testBisection()
{
    double a = 1;
    double b = 2;
    double eps = 1e-8;
    double p;

    p = Bisection(f9, a, b, eps);
    printf("%lf\n", p);
    if(fabs(p - 1.32471795724475) > eps)
    {
        return FAILED;
    }
    return PASSED;
}

static double f10(double x)
{
    return (x + 2) * x * x - 4;
}

static double g10(double x)
{
    return x - f10(x) / (3 * x + 4) / x;
}

/**
 * @brief 
 * @returns 
 * 
 * 
 */
int testPicardRecurtion()
{
    double x = 1;
    double eps = 1e-4;
    double p;

    p = PicardIteration(g10, x, eps);
    printf("%lf\n", p);
    if(fabs(p - 1.130395435) > eps)
    {
        return FAILED;
    }
    return PASSED;
}

static double g11(double x)
{
    return sqrt(4 / (2 + x));
}

/**
 * @brief 
 * @returns 
 * 
 * 
 */
int testSteffensen()
{
    double x = 1.5;
    double eps = 1e-8;
    double p;

    p = SteffensenIteration(g11, x, eps);
    printf("%lf\n", p);
    if(fabs(p - 1.130395435) > eps)
    {
        return FAILED;
    }
    return PASSED;
}

static double f12(double x)
{
    return ((x + 2) * x + 10) * x - 20;
}

static double df12(double x)
{
    return (3 * x + 4) * x + 10;
}

/**
 * @brief 
 * @returns 
 * 
 * 
 */
int testNewtonMethod()
{
    double x = 1;
    double eps = 1e-8;
    double p;

    p = NewtonMethod(f12, df12, x, eps);
    printf("%lf\n", p);
    if(fabs(p - 1.368808108) > eps)
    {
        return FAILED;
    }
    return PASSED;
}

double f13(double x)
{
    return (2 * x * x - 5) * x - 1;
}

/**
 * @brief 
 * @returns 
 * 
 * 
 */
int testSecent()
{
    double x0 = 2;
    double x1 = 1;
    double eps = 1e-8;
    double p;

    p = SecentMethod(f13, x0, x1, eps);
    printf("%lf\n", p);
    if(fabs(p - 1.67298165) > eps)
    {
        return FAILED;
    }
    return PASSED;
}

/**
 * @brief 
 * @returns 
 * 
 * 
 */
int testMuller()
{
    return PASSED;
}
