/** @file testIntegration.c */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "NR.h"
#include "Test.h"

typedef struct _IntTest{
    double (*integrand)(double);
    double (*integral)(double);
    double a;
    double b;
    double TOL;
    double h;
    double hmin;
    double hmax;
}IntTest;

static double integrand1(double x)
{
    return x * sqrt(1 + x * x);
}

static double integral1(double x)
{
    return  pow(1 + x * x, 1.5) / 3;
}

static double integrand2(double x)
{
    return 1 / (1 + x);
}

static double integral2(double x)
{
    return log(x + 1);
}

static IntTest inttest[] = {
    {integrand1, integral1, 0, 3, 1e-10},
    {integrand2, integral2, 0, 3, 1e-10},
    {NULL}
};

int _testInt(double (*f)(IntTest), IntTest t)
{
    double res = f(t);
    double ans = t.integral(t.b) - t.integral(t.a);
    printf("T=%16.12f\tI=%16.12f\tdelta=%16.12f\n", res, ans, fabs(res - ans));
    return (fabs(res - ans) > t.TOL) ? FAILED : PASSED;
}

double _testAdaptiveSimpson(IntTest t)
{
    return AdaptiveSimpsonInt(t.integrand, t.a, t.b, t.TOL);
}

double _testRomberg(IntTest t)
{
    return RombergInt(t.integrand, t.a, t.b, 100, t.TOL);
}

int testAdaptiveSimpson()
{
    for(int i = 0; inttest[i].integrand; i++)
        if(_testInt(_testAdaptiveSimpson, inttest[i]) == FAILED)
            return FAILED;
    return PASSED;
}

int testRomberg()
{
    for(int i = 0; inttest[i].integrand; i++)
        if(_testInt(_testRomberg, inttest[i]) == FAILED)
            return FAILED;
    return PASSED;
}
