
#include <math.h>
#include <stdio.h>
#include "NR.h"

/**
 * @brief Gaussian elimination method to solve linear equation in the form of Ax=b.
 * @param N The rank of the matrix
 * @param A The coefficient matrix
 * @param b The constant vector
 */
double *GaussEli(int N, double **A, double *b)
{
    double *x = (double *)malloc_s(N * sizeof(double));
    double **l = (double **)malloc_s(N * sizeof(double *));
    for(int i = 0; i < N; i++)
    {
        l[i] = (double *)malloc_s(N * sizeof(double));
    }

    for(int k = 0; k < N - 1; k++)
    {
        for(int i = k + 1; i < N; i++)
        {
            l[i][k] = A[i][k] / A[k][k];
            A[i][k] = 0;
            for(int j = k + 1; j < N; j++)
            {
                A[i][j] -= l[i][k] * A[k][j];
            }
            b[i] -= l[i][k] * b[k];
        }
    }

    for(int k = N - 1; k >= 0; k--)
    {
        x[k] = 0;
        for(int j = k + 1; j < N; j++)
        {
            x[k] += A[k][j] * x[j];
        }
        x[k] = (b[k] - x[k]) / A[k][k];
    }

    for(int i = 0; i < N; i++)
    {
        free(*(l + i));
    }
    free(l);
    return x;
}
