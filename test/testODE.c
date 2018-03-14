/** @file testODE.c */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "NR.h"
#include "Test.h"
#include "constants.h"

typedef struct ODETest{
    double (*f)(double, double);
    double (*y)(double);
    double a;
    double b;
    int    N;
    double TOL;
    double h;
    double hmin;
    double hmax;
}ODETest;

static double f1(double t, double y){ return (y * y + y) / t;}
static double y1(double t){ return 2 * t / (1 - 2 * t);}
static double f2(double t, double y){ return -y + t + 1;}
static double y2(double t){ return exp(-t) + t;}
static double f3(double t, double y){ return - y + t * t + 1;}
static double y3(double t){ return - 2 / exp(t) + 3 - 2 * t + t * t;}

static ODETest odetest[] = {
/*    f,  y,   a,   b,  N,  TOL,   h, hmin, hmax */
    {f1, y1, 1.0, 3.0, 20, 1e-2, 0.1            },
    {f2, y2, 1.0, 2.0, 10, 1e-2, 0.1            },
    {f3, y3, 0.0, 1.0, 10, 1e-2, 0.1, 1e-6, 1e-2},
    {NULL}
};

static int _testODE(double *(*f)(ODETest), ODETest t)
{
    double real_y;
    double *result = f(t);
    for(int i = 0; i <= t.N; i ++)
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

static double* _testClassicRK(ODETest t)
{
    return ClassicRungeKutta(t.f, t.a, t.b, t.y(t.a), t.N);
}

int testClassicRK()
{
    for(int i = 0; odetest[i].f; i++)
        if(_testODE(_testClassicRK, odetest[i]) == FAILED)
            return FAILED;
    return PASSED;
}

static double* _testAdamsPECE(ODETest t)
{
    return AdamsPECE(t.f, t.a, t.b, t.f(t.a, t.y(t.a)), t.y(t.a), t.N);
}

int testAdamsPECE()
{
    for(int i = 0; odetest[i].f; i++)
        if(_testODE(_testAdamsPECE, odetest[i]) == FAILED)
            return FAILED;
    return PASSED;
}

static int _testRKF(ODETest t)
{
    ODEsol sol = RKF78(t.f, t.a, t.b, t.y(t.a), 1e-5, t.hmax, t.hmin);
    int step = ODEsolGetStep(sol);
    double *ts = ODEsolGetT(sol);
    double *ys = ODEsolGetY(sol);
    printf("t\t\ty\t\treal_y\n");
    for(int i = 0; i < step; i ++)
    {
        printf("%lf\t%lf\t%lf\n", ts[i], ys[i], t.y(ts[i]));
        if (fabs(ys[i] - t.y(ts[i])) > 1e-5)
        {
            return FAILED;
        }
    }
    return PASSED;
}

int testRKF()
{
    return _testRKF(odetest[2]);
}

static double *Sf1(double t, double *y)
{
    double *f = newArray1d(2);
    f[0] = -4 * y[0] - 2 * y[1] + cos(t) + 4 * sin(t);
    f[1] = 3 * y[0] + y[1] - 3 * sin(t);
    return f;
}

static double *Sy1(double x)
{
	double *y = newArray1d(2);
    y[0] = exp(-2 * x) * (-2 + 2 * exp(x) + exp(2 * x) * sin(x));
    y[1] = -exp(-2 * x) * (3 * exp(x) - 2);
    return y;
}

typedef struct SODETest{
    double* (*f)(double, double*);
    double* (*y)(double);
    double a;
    double b;
    int    N;
    int    m;
    double TOL;
    double h;
    double hmin;
    double hmax;
}SODETest;

static SODETest sodetest[] = {
/*     f,   y,   a,   b,  N,  m,  TOL,   h, hmin, hmax */
    {Sf1, Sy1, 0.0, 1.0, 10,  2, 1e-8, 0.1, 1e-4, 1.0},
    {Sf1, Sy1, 0.0, 1.0, 10,  2, NAN , 0.1, 1e-4, 1.0},
    {NULL}
};

/**
 * @brief
 * @returns
 */
int _testSODE(SODEsol (*f)(SODETest t), SODETest t)
{
    SODEsol sol = f(t);

    double *ts = SODEsolGetT(sol);
    double **ys = SODEsolGetY(sol);
    int steps = SODEsolGetStep(sol);

    printf("%8s%16s%16s%16s%16s\n", "t", "result1", "y1", "result2", "y2");
    for(int i = 1; i < steps; i ++)
    {
		double *ans = t.y(ts[i]);
        printf("%8.4lf", ts[i]);
        for(int j = 0 ; j < t.m; j++)
        {
            printf("%16.11f%16.11f", ys[i][j], ans[j]);
            if(!isnan(t.TOL) && 
             && fabs(ys[i][j] - ans[j])/(ts[i] - ts[i - 1]) > t.TOL * i)
            {
                return FAILED;
            }
        }
        printf("\n");
        delArray1d(ans);
    }
    delArray1d(ts);
    delArray2d(ys, steps);
    return PASSED;
}

SODEsol _testSODERKF78(SODETest t)
{
    return SODERKF78(t.f, t.y(t.a), t.a, t.b, t.m, t.h, t.TOL, t.hmax, t.hmin);
}

SODEsol _testSODERKF45(SODETest t)
{
    return SODERKF45(t.f, t.y(t.a), t.a, t.b, t.m, t.h, t.TOL, t.hmax, t.hmin);
}

SODEsol _testSODERungeKutta(SODETest t)
{
    return SODERungeKutta(t.f, t.a, t.b, t.y(t.a), t.m, t.N);
}

int testSODERKF78()
{
    return _testSODE(_testSODERKF78, sodetest[0]);
}

int testSODERKF45()
{
    return _testSODE(_testSODERKF45, sodetest[0]);
}

int testSODERungeKutta()
{
    return _testSODE(_testSODERungeKutta, sodetest[1]);
}
