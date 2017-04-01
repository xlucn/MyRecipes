/**
 * @file 
 * @brief test solving ODEs
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "NR.h"
#include "Test.h"

typedef struct _ODETest{
    double (*f)(double, double);
    double (*y)(double);
    double a;
    double b;
    double TOL;
    double h;
    double hmin;
    double hmax;
}ODETest;

typedef struct _SODETest{
    double* (*f)(double, double);
    double (*y)(double);
    double a;
    double b;
    double TOL;
    double h;
    double hmin;
    double hmax;
}SODETest;

static double f1(double t, double y)
{
    return (y * y + y) / t;
}

static double y1(double t)
{
    return 2 * t / (1 - 2 * t);
}

static double f2(double t, double y)
{
    return -y + t + 1;
}

static double y2(double t)
{
    return exp(-t) + t;
}

static double f3(double t, double y)
{
    return - y + t * t + 1;
}

static double y3(double t)
{
    return - 2 / exp(t) + 3 - 2 * t + t * t;
}

static ODETest odetest[] = {
    {f1, y1, 1.0, 3.0, 0.01, 0.5},
    {f2, y2, 1.0, 2.0, 0.01, 0.2},
    {f3, y3, 0.0, 1.0, 0.01, 0.2},
    {NULL}
};

static int _testClassicRK(ODETest t)
{
    printf("Testing Classic Runge-Kutta method\n");
    double real_y;
    int N = (int)((t.b - t.a) / t.h + 0.5);
    double *result = ClassicRungeKutta(t.f, t.a, t.b, t.y(t.a), N);
    for(int i = 0; i < N + 1; i ++)
    {
        real_y = t.y(t.a + i * t.h);
        printf("%lf\t%lf\n", result[i], real_y);
        if (fabs(result[i] - real_y) > t.TOL)
        {
            return FAILED;
        }
    }
    return PASSED;
}

int testClassicRK()
{
    for(int i = 0; odetest[i].f; i++)
        if(_testClassicRK(odetest[i]) == FAILED)
            return FAILED;
    return PASSED;
}

static int _testAdamsPECE(double a, double b, int N, double f(double, double), double y(double), double eps)
{
    double y0 = y(a);
    double dy0 = f(a, y0);
    double *ans = (double*)malloc_s((N + 1) * sizeof(double));
    for(int i = 0; i < N + 1; i++)
    {
        ans[i] = y(a + (b - a) * i / N);
    }

    double *result = AdamsPECE(f, a, b, dy0, y0, N);
    printf("Testing Adams PECE\n");
    for(int i = 0; i < N + 1; i ++)
    {
        printf("%lf\t%lf\n", result[i], ans[i]);
        if(fabs(result[i] - ans[i]) > eps)
        {
            return FAILED;
        }
    }
    return PASSED;
}

/**
 * @brief 
 * @returns 
 */
int testAdamsPECE()
{
    double eps = 1e-3;
    return _testAdamsPECE(0, 1, 10, f2, y2, eps)
        || _testAdamsPECE(1, 2, 20, f1, y1, eps)
        || _testAdamsPECE(0, 1, 10, f3, y3, eps);
}

static int _testRKF(double a, double b, double TOL, double hmax, double hmin, double f(double, double), double y(double))
{
    double y0 = y(a);
    double *result;

    printf("Testing RKF Method\n");
    result = RKF78(f, a, b, y0, TOL, hmax, hmin);
    if(result == NULL)
    {
        fprintf(stderr, "RKF method got NULL result.\n");
        return FAILED;
    }
    printf("t\t\th\t\ty\t\treal_y\n");
    for(int i = 0; i < result[0]; i ++)
    {
        for(int j = 0; j < 3; j++)
        {
            printf("%lf\t", result[i * 3 + j + 1]);
        }
        printf("%lf\n", y(result[i * 3 + 1]));

        // Check the results
        if (fabs(result[i * 3 + 3] - y(result[i * 3 + 1])) > 1e-5)
        {
            return FAILED;
        }
    }
    return PASSED;
}

/**
 * @brief 
 * @returns 
 * 
 * 
 */
int testRKF()
{
    return _testRKF(0, 1, 1e-5, 1e-2, 1e-6, f3, y3);
}

static double *Sf1(double t, double *y)
{
    static double f[2];
    f[0] = -4 * y[0] - 2 * y[1] + cos(t) + 4 * sin(t);
    f[1] = 3 * y[0] + y[1] - 3 * sin(t);
    return f;
}

static double *Sy1(double x)
{
	static double y[2];
    y[0] = exp(-2 * x) * (-2 + 2 * exp(x) + exp(2 * x) * sin(x));
    y[1] = -exp(-2 * x) * (3 * exp(x) - 2);
    return y;
}

/**
 * @brief 
 * @returns 
 */
int testSODERungeKutta()
{
    double *(*f)(double,double*) = Sf1;
    double *(*y)(double) = Sy1;
    double a = 0;
    double b = 1;
    double y0[] = {0, -1};
    //int N = 10;
    int m = 2;
    SODEsol sol = SODERKF(f, y0, a, b, m, 0.1, 1e-8, 1e4, 1e-4, 13);
    //SODERungeKutta(f, a, b, y0, m, N);

    double *t = SODEsolGetT(sol);
    double **res = SODEsolGetY(sol);
    int steps = SODEsolGetStep(sol);


    printf("Testing Runge-Kutta Method for a System of ODEs\n");
    printf("%8s%16s%16s%16s%16s\n", "t", "result1", "y1", "result2", "y2");
    for(int i = 0; i < steps; i ++)
    {
		double *ans = y(t[i]);
        printf("%8.2lf", t[i]);
        for(int j = 0 ; j < m; j++)
        {
            printf("%16lf%16lf", res[i][j], ans[j]);
        }
        printf("\n");
    }
    free(t);
    for (int i = 0; i < steps; i++)
	{
		free(res[i]);
	}
	free(res);
	
    return PASSED;
}

