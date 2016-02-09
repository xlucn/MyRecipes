//Author LuXu
//Basic functions used by more than one task

#include <math.h>
#include <stdlib.h>
#include <NumericalRecipes.h>
#include <LibFunction.h>

double DividedDiff(double f(double), double *x, int N)
{
    return N == 1 ? f(*x) : (DividedDiff(f, x + 1, N - 1) - DividedDiff(f, x, N - 1)) / (x[N - 1] - x[0]);
}

/*
All the Divided differences of a array up to order k
*/
double **FullDividedDiff(double f(double), double *x, int N, int k)
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
/*

*/
double **DividedDiffMatrix(double f(double), double *x, int N)
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

/*
拉格朗日多项式
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
/*
Chebyshev ploy
*/
double Chebyshev(int n, double x)
{
    return cos(n * acos(x));
}
