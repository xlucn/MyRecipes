//LuXu
//Caculating integrals

#include <stdlib.h>
#include <math.h>
#include "NR.h"
#include "LibFunction.h"

double TrapezoidalInt(double f(double), double a, double b)
{
    return (b - a) / 2 * (f(a) + f(b));
}

double SimpsonInt(double f(double), double a, double b)
{
    return (b - a) / 6 * (f(a) + 4 * f((a + b) / 2) + f(b));
}

double CompositeSimpsonInt(double f(double), double a, double b, int N)
{
    double Integral = 0;
    double h = (b - a) / N / 2;

    for (int i = 0; i < N; i++)
    {
        Integral += f(a + 2 * i * h) * h / 3;
        Integral += f(a + (2 * i + 1) * h) * 4 * h / 3;
        Integral += f(a + (2 * i + 2) * h) * h / 3;
    }

    return Integral;
}

double RombergInt(double f(double), double a, double b, int N, double eps)
{
    double result;
    double h = b - a;
    double **T = (double**)malloc_s(N * sizeof(double*));

    for (int i = 0; i < N; i++)
    {
        *(T + i) = (double*)malloc_s((i + 1) * sizeof(double));
        if (i == 0)
        {
            //T[0][0]
            T[0][0] = (f(b) + f(a)) * h / 2;
        }
        else
        {
            //T[i][0]
            double temp = 0;
            for (int m = 0; m < pow(2, i - 1); m++)
            {
                temp += f(a + (m + 0.5) * h);
            }
            T[i][0] = (T[i - 1][0] + h * temp) / 2;

            //T[i][j]
            for (int j = 1; j < i + 1; j++)
            {
                T[i][j] = (pow(4, j) * T[i][j - 1] - T[i - 1][j - 1]) / (pow(4, j) - 1);
            }

            if (fabs(T[i][i] - T[i][i - 1]) <= eps)
            {
                return T[i][i];
            }

            free(*(T + i - 1));

            h = h / 2;
        }
    }
    result = T[N][N - 1];
    free(*(T + N - 1));
    free(*T);
    return result;

}

/*
自适应Simpson积分，利用迭代法
*/
double AdaptiveSimpsonInt(double f(double), double a, double b, double TOL)
{
    double Stotal = SimpsonInt(f, a, b);
    double Sleft  = SimpsonInt(f, a, (a + b) / 2);
    double Sright = SimpsonInt(f, (a + b) / 2, b);

    if (fabs(Stotal - Sleft - Sright) > 15 * TOL)
    {
        Sleft = AdaptiveSimpsonInt(f, a, (a + b) / 2, TOL / 2);
        Sright = AdaptiveSimpsonInt(f, (a + b) / 2, b, TOL / 2);
    }

    return Sleft + Sright;
}
