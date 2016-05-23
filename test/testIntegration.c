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
	{integrand1, integral1, 0, 3, 1e-4},
	{integrand2, integral2, 0, 3, 1e-4}
};

int _testAdaptiveSimpson(IntTest item)
{
    double res = AdaptiveSimpsonInt(item.integrand, item.a, item.b, item.TOL);
	double ans = item.integral(item.b) - item.integral(item.a);
	printf("T=%.12f\nI=%.12f\ndelta=%.12f\n", res, ans, fabs(res - ans));
	return (fabs(res - ans) > item.TOL) ? FAILED : PASSED;
}

int _testRomberg(IntTest item)
{
    int N = 100;
    double TOL = 1e-10;
	double res = RombergInt(item.integrand, item.a, item.b, N, TOL);
	double ans = item.integral(item.b) - item.integral(item.a);
	printf("T=%.12f\nI=%.12f\ndelta=%.12f\n", res, ans, fabs(res - ans));
	return (fabs(res - ans) > TOL) ? FAILED : PASSED;
}


/**
 * @brief test adaptive Simpson method
 * @returns PASSED or FAILED
 */
int testAdaptiveSimpson()
{
	return _testAdaptiveSimpson(inttest[0])
		|| _testAdaptiveSimpson(inttest[1]);
}


/**
 * @brief test Romberg method
 * @returns PASSED or FAILED
 */
int testRomberg()
{
	return _testRomberg(inttest[0])
		|| _testRomberg(inttest[1]);
}
