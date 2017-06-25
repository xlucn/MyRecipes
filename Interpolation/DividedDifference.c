/** @file DividedDifference.c */
#include "NR.h"
/**
 * @brief Divided difference of function f on nodes x[n]
 * @param f function
 * @param x nodes
 * @param N number of nodes
 * @returns divided difference of function f on nodes x[n]
 *
 * all orders of divided differences are defined as follows:
 * f[x] = f(x),
 * f[x1, x2] = (f(x2)-f(x1))/(x2-x1),
 * ......
 * f[x1, x2, ... , xn] = (f[x2, x3, ..., xn]-f[x1, x2, ..., x(n-1)])/(xn-x1)
 */
double DividedDiff(double (*f)(double), double *x, int N)
{
    double d = 0;
    for(int i = 0; i < N; i++)
    {
        double fi = f(x[i]);
        for(int j = 0; j < N; j++) if(j != i)
            fi /= (x[i] - x[j]);
        d += fi;
    }
    return d;
}

/**
 * @brief All the Divided differences of a array up to order k
 * @param f function
 * @param x nodes
 * @param N number of nodes
 * @param k upper order of divided difference
 * @returns All the Divided differences of a array up to order k
 */
double **FullDividedDiff(double (*f)(double), double *x, int N, int k)
{
    double **d = newArray2d(k, N);
    for (int j = 0; j < k; j++)
    {
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
 * @brief calculate a matrix of the divided differnces.
 * @param f function
 * @param x nodes
 * @param N number of nodes
 * @returns a matrix of the divided differnces.
 *
 * d[i][j] is f[xi, ..., xj] as in DividedDiff()
 */
double **DividedDiffMatrix(double (*f)(double), double *x, int N)
{
    double **d = newArray2d(N, N);

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

