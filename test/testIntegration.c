/**
 * @file testIntegration.c
 * @brief testIntegration.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "NR.h"
#include "Test.h"

/**
 * @brief Integral test unit
 */
typedef struct _IntTest{
    double (*integrand)(double); /**< empty */
    double (*integral)(double); /**< empty */
    double a; /**< empty */
    double b; /**< empty */
    double TOL; /**< empty */
    double h; /**< empty */
    double hmin; /**< empty */
    double hmax; /**< empty */
}IntTest;

static double integrand1(double x) { return x * sqrt(1 + x * x); }
static double integral1(double x)  { return  pow(1 + x * x, 1.5) / 3; }
static double integrand2(double x) { return 1 / (1 + x); }
static double integral2(double x)  { return log(x + 1); }

static IntTest inttest[] = {
    {integrand1, integral1, 0, 3, 1e-10},
    {integrand2, integral2, 0, 3, 1e-10},
    {NULL}
};

static int _testInt(double (*f)(IntTest), IntTest t)
{
    double res = f(t);
    double ans = t.integral(t.b) - t.integral(t.a);
    printf("T=%16.12f\tI=%16.12f\tdelta=%16.12f\n", res, ans, fabs(res - ans));
    return (fabs(res - ans) > t.TOL) ? FAILED : PASSED;
}

static double _testAdaptiveSimpson(IntTest t)
{
    return AdaptiveSimpsonInt(t.integrand, t.a, t.b, t.TOL);
}

static double _testRomberg(IntTest t)
{
    return RombergInt(t.integrand, t.a, t.b, 100, t.TOL);
}

/**
 * @brief testAdaptiveSimpson
 * @return integer if test passed
 */
int testAdaptiveSimpson()
{
    for(int i = 0; inttest[i].integrand; i++)
        if(_testInt(_testAdaptiveSimpson, inttest[i]) == FAILED)
            return FAILED;
    return PASSED;
}

/**
 * @brief testRomberg
 * @return integer if test passed
 */
int testRomberg()
{
    for(int i = 0; inttest[i].integrand; i++)
        if(_testInt(_testRomberg, inttest[i]) == FAILED)
            return FAILED;
    return PASSED;
}
