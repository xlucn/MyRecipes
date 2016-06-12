#include <math.h>
#include <stdio.h>
#include "NR.h"

/**
 * @brief Gaussian elimination method to solve linear equation in the form of Ax=b.
 * @param N The numbers of the variables
 * @param a The coefficient matrix in a 1-dimension array
 * @param b The constant vector
 * @return An array of numbers, NULL if the equation has no solution
 */
double *GaussEli(int N, double *a, double *b)
{
    double **A = (double**)malloc_s(N * sizeof(double*));
    for (int i = 0; i < N; i++)
    {
        A[i] = (double*)malloc_s((N + 1) * sizeof(double));
        for (int j = 0; j < N; j++)
        {
            A[i][j] = a[i * N + j];
        }
        A[i][N] = b[i];
    }
    
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
            for(int j = k + 1; j < N + 1; j++)
            {
                A[i][j] -= l[i][k] * A[k][j];
            }
        }
    }

    for(int k = N - 1; k >= 0; k--)
    {
        double t = 0;
        for(int j = k + 1; j < N; j++)
        {
            t += A[k][j] * x[j];
        }
        x[k] = (A[k][N] - t) / A[k][k];
    }

    for(int i = 0; i < N; i++)
    {
        free(*(l + i));
    }
    free(l);
    return x;
}
