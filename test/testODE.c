#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "NR.h"
#include "Test.h"

/**
 * testing solving ODEs
 */

static double ODEf1(double t, double y)
{
    return (y * y + y) / t;
}

static double ODEy1(double t)
{
    return 2 * t / (1 - 2 * t);
}

static double ODEf2(double t, double y)
{
    return -y + t + 1;
}

static double ODEy2(double t)
{
    return exp(-t) + t;
}

static double ODEf3(double t, double y)
{
    return - y + t * t + 1;
}

static double ODEy3(double t)
{
    return - 2 / exp(t) + 3 - 2 * t + t * t;
}

static int _testClassicRK(double a, double b, int N, double ODEf(double, double), double ODEy(double), double eps)
{
    double t, real_y;
    double y0 = ODEy(a);

    double *result = ClassicRungeKutta(ODEf, a, b, y0, N);
    printf("Testing Classic Runge-Kutta method\n");
    for(int i = 0; i < N + 1; i ++)
    {
        t = a + i * (b - a) / N;
        real_y = ODEy(t);
        printf("%lf\t%lf\n", result[i], real_y);
        if (fabs(result[i] - real_y) > eps)
        {
            fprintf(stderr, "Classic Runge-Kutta method accuracy not enough\n");
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
int testClassicRK()
{
    double eps = 0.01;
    return _testClassicRK(1, 3, 4, ODEf1, ODEy1, eps)
        || _testClassicRK(1, 2, 5, ODEf2, ODEy2, eps)
        || _testClassicRK(0, 1, 5, ODEf3, ODEy3, eps);
}


/**
 * @brief 
 * @param a 
 * @param b 
 * @param N 
 * @param ODEf 
 * @param ODEy 
 * @param eps 
 * @returns 
 * 
 * 
 */
static int _testAdamsPECE(double a, double b, int N, double ODEf(double, double), double ODEy(double), double eps)
{
    double y0 = ODEy(a);
    double dy0 = ODEf(a, y0);
    double *ans = (double*)malloc_s((N + 1) * sizeof(double));
    for(int i = 0; i < N + 1; i++)
    {
        ans[i] = ODEy(a + (b - a) * i / N);
    }

    double *result = AdamsPECE(ODEf, a, b, dy0, y0, N);
    printf("Testing Adams PECE\n");
    for(int i = 0; i < N + 1; i ++)
    {
        printf("%lf\t%lf\n", result[i], ans[i]);
        if(fabs(result[i] - ans[i]) > eps)
        {
            fprintf(stderr, "Adams PECE method accuracy not enough\n");
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
int testAdamsPECE()
{
    double eps = 1e-3;
    return _testAdamsPECE(0, 1, 10, ODEf2, ODEy2, eps)
        || _testAdamsPECE(1, 2, 20, ODEf1, ODEy1, eps)
        || _testAdamsPECE(0, 1, 10, ODEf3, ODEy3, eps);
}

/**
 * @brief 
 * @param a 
 * @param b 
 * @param TOL 
 * @param hmax 
 * @param hmin 
 * @param ODEf 
 * @param ODEy 
 * @returns 
 * 
 * 
 */
static int _testRKF(double a, double b, double TOL, double hmax, double hmin, double ODEf(double, double), double ODEy(double))
{
    double y0 = ODEy(a);
    double *result;

    printf("Testing RKF Method\n");
    result = RKF78(ODEf, a, b, y0, TOL, hmax, hmin);
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
        printf("%lf\n", ODEy(result[i * 3 + 1]));

        // Check the results
        if (fabs(result[i * 3 + 3] - ODEy(result[i * 3 + 1])) > 1e-5)
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
    return _testRKF(0, 1, 1e-5, 1e-2, 1e-6, ODEf3, ODEy3);
}

static double *SODEf1(double t, double *y)
{
    double *f = (double*)malloc_s(2 * sizeof(double));
    f[0] = -4 * y[0] - 2 * y[1] + cos(t) + 4 * sin(t);
    f[1] = 3 * y[0] + y[1] - 3 * sin(t);
    return f;
}

static double *SODEy1(double x)
{
	double *y = (double*)malloc_s(2 * sizeof(double));
    y[0] = exp(-2 * x) * (-2 + 2 * exp(x) + exp(2 * x) * sin(x));
    y[1] = -exp(-2 * x) * (3 * exp(x) - 2);
    return y;
}

/**
 * @brief 
 * @returns 
 * 
 * 
 */
int testSODERungeKutta()
{
    double *(*f)(double,double*) = SODEf1;
    double *(*y)(double) = SODEy1;
    double a = 0;
    double b = 1;
    double y0[] = {0, -1};
    //int N = 10;
    int m = 2;
    SODEsol sol = SODERKF(f, y0, a, b, m, 0.1, 1e-8, 1e4, 1e-4, 13);
    //SODERungeKutta(f, a, b, y0, m, N);

    double *t = sol.t;
    double **res = sol.y;
    int steps = sol.step;
    

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

