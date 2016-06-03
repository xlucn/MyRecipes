//Author LuXu
//Basic functions used by more than one task

#include <stdlib.h>
#include <math.h>
#include "NR.h"
/**
 * @brief Divided difference of function f on nodes x[n]
 * f[x] = f(x),
 * f[x1, x2] = (f(x2)-f(x1))/(x2-x1),
 * ......
 * f[x1, x2, ... , xn] = (f[x2, x3, ..., xn]-f[x1, x2, ..., x(n-1)])/(xn-x1)
 */
double DividedDiff(double(*f)(double), double *x, int N)
{
    return N == 1 ? f(*x) : (DividedDiff(f, x + 1, N - 1) - DividedDiff(f, x, N - 1)) / (x[N - 1] - x[0]);
}
/**
 * @brief All the Divided differences of a array up to order k
 */
double **FullDividedDiff(double(*f)(double), double *x, int N, int k)
{
    double **d = (double**)malloc_s(k * sizeof(double*));
    for (int j = 0; j < k; j++)
    {
        d[j] = (double*)malloc_s((N - j) * sizeof(double));
        for (int i = 0; i < N - j; i++)
        {
            if (j == 0)
            {
                d[0][i] = f(x[i]);
            }
            else
            {
                d[j][i] = (d[j - 1][i + 1] - d[j - 1][i]) / (x[i + j] - x[i]);
            }
        }
    }
    return d;
}
/**
 * @brief return a matrix of the divided differnces. d[i][j] is f[xi, ..., xj]
 */
double **DividedDiffMatrix(double(*f)(double), double *x, int N)
{
    double **d = (double**)malloc_s(N * sizeof(double*));
    for (int i = 0; i < N; i++)
    {
        d[i] = (double*)malloc_s(N * sizeof(double));
    }

    for (int j = 0; j < N; j++)
    {
        for (int i = N - 1; i >= 0; i--)
        {
            if (i > j)
            {
                d[i][j] = 0;
            }
            else if (i == j)
            {
                d[i][j] = f(x[i]);
            }
            else
            {
                d[i][j] = (d[i + 1][j] - d[i][j - 1]) / (x[j] - x[i]);
            }
        }
    }
    return d;
}

/**
 * @brief Lagrange polynomial
 * 输入N和插值点向量a[N]，以及变量值x，返回拉格朗日基本多项式的值
 */
double *LagrangePoly(double *a, double x, int N)
{
    double *l = (double *)malloc_s(N * sizeof(double));
    for(int i = 0; i < N; i++)
    {
        l[i] = 1;
        for (int j = 0; j < N; j++)
        {
            if (j != i)
            {
                l[i] = l[i] * (x - a[j]) / (a[i] - a[j]);
            }
        }
    }
    return l;
}
/**
 * @brief Chebyshev polynomial(first kind)
 */
double Chebyshev(int n, double x)
{
    return cos(n * acos(x));
}
