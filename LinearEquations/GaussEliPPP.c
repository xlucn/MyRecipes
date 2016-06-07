
#include <math.h>
#include <stdio.h>
#include "NR.h"

/**
 * @brief Gauss elimination with partial pivoting proportionally
 * @param N the rank of the equation
 * @param a the augmented matrix
 */
double *GaussEliPPP(int N, double **a)
{
    int r = 0;
    double q, *t; // for swap
    double *max = (double*)malloc_s(N * sizeof(double));
    double *x = (double*)malloc_s(N * sizeof(double));

    // get the largest number of each line
    for (int i = 0; i < N; i++)
    {
        max[i] = 0;
        for (int j = 0; j < N; j++)
        {
            max[i] = (fabs(max[i]) < fabs(a[i][j])) ? fabs(a[i][j]) : max[i];
        }
        if (max[i] == 0)
        {
            printf("A is singular\n");
        }
    }

    // do the gauss elimination with partial pivoting proportionally
    for (int k = 0; k < N - 1; k++)
    {
        // choose the pivot element
        for (int i = k; i < N; i++)
        {
            r = (fabs(a[r][k] / max[r]) < fabs(a[i][k] / max[i])) ? i : r;
        }
        if (a[r][k] == 0)
        {
            printf("A is singular\n");
            return NULL;
        }
        if (r != k)
        {
            q = max[k]; max[k] = max[r]; max[r] = q;
            t = a[k]; a[k] = a[r]; a[r] = t;
        }

        for (int i = k + 1; i < N; i++)
        {
            a[i][k] = a[i][k] / a[k][k];
            for (int j = k + 1; j < N + 1; j++)
            {
                a[i][j] -= a[i][k] * a[k][j];
            }
        }
    }

    if (a[N - 1][N - 1] == 0)
    {
        printf("A is singular\n");
        return NULL;
    }

    for (int k = N - 1; k >= 0; k--)
    {
        x[k] = 0;
        for (int j = k + 1; j < N; j++)
        {
            x[k] += a[k][j] * x[j];
        }
        x[k] = (a[k][N] - x[k]) / a[k][k];
    }
    free(max);
    return x;
}

/**
 * @brief Gauss Jordan elimination method to solve a system of linear equations
 * @param N the rank of the equation
 * @param a the augmented matrix
 */
double *GaussJordanEli(int N, double **a)
{
    int* max = (int*)malloc_s(N * sizeof(int));
    double* x = (double*)malloc_s(N * sizeof(double));

    for (int k = 0; k < N; k++)
    {
        max[k] = 0;
        for (int i = 0; i < N; i++)
        {
            int flag = 1;
            for (int j = 0; j < k; j++)
            {
                if (i == max[j])
                {
                    flag = 0;
                }
            }
            if (flag)
            {
                max[k] = (fabs(a[max[k]][k]) < fabs(a[i][k])) ? i : max[k];
            }
        }
        if (a[max[k]][k] == 0)
        {
            printf("A is singular\n");
            return NULL;
        }

        for (int i = 0; i < N; i++)
        {
            if (i != max[k])
            {
                a[i][k] = a[i][k] / a[max[k]][k];
                for (int j = k + 1; j < N + 1; j++)
                {
                    a[i][j] -= a[i][k] * a[max[k]][j];
                }
            }

        }
    }

    for (int k = 0; k < N; k++)
    {
        x[k] = a[max[k]][N] / a[max[k]][k];
    }
    free(max);
    return x;
}
