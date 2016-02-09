//2015.11.27
//卢旭
//Solving Ordinary differential equations
#include <stdlib.h>
#include <math.h>
#include <NumericalRecipes.h>
#include <LibFunction.h>

/*
Euler Method to solve ODE. y0:initial value, f(t, y) = dy/dt
*/
double* Euler(double f(double, double), double a, double b, double y0, int N)
{
    double h = (b - a) / N;
    double x = a;
    double* y = (double*)malloc_s((N + 1) * sizeof(double));
    y[0] = y0;

    for(int i = 0; i < N; i++)
    {
        y[i + 1] += h * f(x, y[i]);
        x += h;
    }
    return y;
}

/*
2阶方法的一般步骤
*/
static double* TwoStageRungeKutta(int N, double y0, double a, double b, double f(double, double), double para)
{
    double c1 = 1 -  0.5 / para;
    double c2 = 1 - c1;
    double k1;
    double k2;
    double h = (b - a) / N;
    double x = a;
    double y = y0;
    double* result = (double*)malloc_s(N * sizeof(double));

    for(int i = 0; i < N; i++)
    {
        k1 = f(x, y);
        k2 = f(x + para * h, y + para * h * k1);
        y += h * (c1 * k1 + c2 * k2);
        x += h;
        result[i] = y;
    }
    return result;
}
/*
改善的欧拉方法
*/
double* ImprovedEuler(double f(double, double), double a, double b, double y0, int N)
{
    return TwoStageRungeKutta(N, y0, a, b, f, 1);
}
/*
中点方法或变形的欧拉方法
*/
double* MID(double f(double, double), double a, double b, double y0, int N)
{
    return TwoStageRungeKutta(N, y0, a, b, f, 1 / 2);
}
/*
Heun方法
*/
double* Heun(double f(double, double), double a, double b, double y0, int N)
{
    return TwoStageRungeKutta(N, y0, a, b, f, 2 / 3);
}

/*
三阶Heun方法
*/
double *ThreeStageHeun(int N, double y0, double a, double b, double f(double, double))
{
    return 0;
}

/*
ThreeStageRungeKuttaMathod
*/
double *ThreeStageRungeKuttaMathod(double f(double, double), double a, double b, double y0, int N)
{
    return 0;
}

/*
Classic Runge-Kutta Method
*/
double *ClassicRungeKutta(double f(double, double), double a, double b, double y0, int N)
{
    double k[4];
    double h = (b - a) / N;
    double x = a;
    double y = y0;
    double* result = (double*)malloc_s((N + 1) * sizeof(double));
    result[0] = y0;

    for(int i = 0; i < N; i++)
    {
        k[0] = f(x, y);
        k[1] = f(x + h / 2, y + h * k[0] / 2);
        k[2] = f(x + h / 2, y + h * k[1] / 2);
        k[3] = f(x + h, y + h * k[2]);
        y += h / 6 * (k[0] + 2 * k[1] + 2 * k[2] + k[3]);
        x += h;
        result[i + 1] = y;
    }

    return result;
}

/*
Runge-Kutta-Fehlberg Method, one of the adaptive Runge-Kutta methods.
*/
double *RKF(double f(double,double), double a, double b, double y0, double TOL, double hmax, double hmin)
{
    static double A[6][6] =
    {
        {0},
        {1.0/4.0, 0},
        {3.0/32.0, 9.0/32.0, 0},
        {1932.0/2197.0, -7200.0/2197.0, 7296.0/2197.0, 0},
        {439.0/216.0, -8.0, 3680.0/513.0, -845.0/4104.0, 0},
        {-8.0/27.0, 2.0, -3544.0/2565.0, 1859.0/4104.0, -11.0/40.0, 0}
    };
    static double B[6] = {25.0/216.0, 0, 1408.0/2565.0, 2197.0/4104.0, -1.0/5.0, 0};
    static double Bstar[6] = {16.0/135.0, 0, 6656.0/12825.0, 28561.0/56430.0, -9.0/50.0, 2.0/55.0};
    static double C[6] = {0, 1.0/4.0, 3.0/8.0, 12.0/13.0, 1, 1.0/2.0};

    int step = 0; //the total steps
    double t = a;
    double y = y0;
    double h = hmax;
    double k[6];
    double R = 0;
    double delta;
    double *result;
    step++;
    result = (double*)malloc_s((step * 3 + 1) * sizeof(double));
    result[1] = t;
    result[2] = h;
    result[3] = y;

    while (t < b)
    {
        for(int i = 0; i < 6; i++)
        {
            y = result[3 * step];
            for(int j = 0; j < i; j++)
            {
                y += A[i][j] * k[j];
            }
            k[i] = h * f(t + C[i] * h, y);
            R += (B[i] - Bstar[i]) * k[i] / h;
        }
        y = result[3 * step];
        R = fabs(R);
        delta = pow((TOL / R / 2.0), 0.25);

        if (R <= TOL)
        {
            t += h;
            for(int i = 0; i < 5; i++)
            {
                y += B[i] * k[i];
            }
            step++;
            result = (double*)realloc(result, (step * 3 + 1) * sizeof(double));
            result[step * 3 - 2] = t;
            result[step * 3 - 1] = h;
            result[step * 3] = y;
        }

        h = delta < 0.1 ? 0.1 * h : (delta > 4 ? 4 * h: delta * h);

        if (h >= hmax)
        {
            h = hmax;
        }
        if (h < hmin)
        {
            return NULL;
        }
    }
    result[0] = step;
    return result;
}

/*
Adams显式和隐式方法的PECE模式校正方法，这里k=1，用经典Runge-Kutta方法提供初值
*/
double *AdamsPECE(double f(double, double), double a, double b, double dy0, double y0, int N)
{
    double *y = (double *)malloc_s((N + 1) * sizeof(double));
    double *dy = (double *)malloc_s((N + 1) * sizeof(double));
    double h = (b - a) / N;
    y[0] = y0;
    y[1] = ClassicRungeKutta(f, a, b, y0, N)[1];
    dy[0] = dy0;
    dy[1] = f(a + h, y[1]);

    for(int i = 1; i < N; i ++)
    {
        y[i + 1] = y[i] + h / 2 * (3 * dy[i] - dy[i - 1]);
        dy[i + 1] = f(a + (i + 1) * h, y[i + 1]);
        y[i + 1] = y[i] + h / 2 * (dy[i] + dy[i + 1]);
        dy[i + 1] = f(a + (i + 1) * h, y[i + 1]);
    }
    free(dy);
    return y;
}

double **SODERungeKutta(double (**f)(double, double*), double a, double b, double *y0, int m, int N)
{
    double h = (b - a) / N;
    double t = a;
    double *w = (double *)malloc_s(m * sizeof(double));
    double **y = (double**)malloc_s((N + 1) * sizeof(double*));
    for(int i = 0; i < N + 1; i++)
    {
        *(y + i) = (double*)malloc_s(m * sizeof(double));
    }
    double **k = (double**)malloc_s(4 * sizeof(double*));
    for(int i = 0; i < 4; i++)
    {
        *(k + i) = (double*)malloc_s(m * sizeof(double));
    }

    for(int i = 0; i < m; i++)
    {
        w[i] = y0[i];
        y[0][i] = w[i];
    }
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < m; j++)
        {
            k[0][j] = h * f[j](t, w);
            w[j] = y[i][j] + k[0][j] / 2;
        }
        for(int j = 0; j < m; j++)
        {
            k[1][j] = h * f[j](t + h / 2, w);
            w[j] = y[i][j] + k[1][j] / 2;
        }
        for(int j = 0; j < m; j++)
        {
            k[2][j] = h * f[j](t + h / 2, w);
            w[j] = y[i][j] + k[2][j];
        }
        for(int j = 0; j < m; j++)
        {
            k[3][j] = h * f[j](t + h, w);
            y[i + 1][j] = y[i][j] + (k[0][j] + k[1][j] * 2 + k[2][j] * 2 + k[3][j]) / 6;
        }
        t += h;
    }
    return y;
}
