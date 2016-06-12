#include <math.h>
#include <stdio.h>
#include "NR.h"

/**
 * @brief Gauss elimination with partial pivoting proportionally
 * @param N The numbers of the variables
 * @param a The coefficient matrix in a 1-dimension array
 * @param b The constant vector
 * @return An array of numbers, NULL if the equation has no solution
 */
double *GaussEliPPP(int N, double *a, double *b)
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
    
    double *max = (double*)malloc_s(N * sizeof(double));
    double *x = (double*)malloc_s(N * sizeof(double));

    /* get the largest number of each line */
    for (int i = 0; i < N; i++)
    {
        max[i] = 0;
        for (int j = 0; j < N; j++)
        {
            max[i] = (fabs(max[i]) < fabs(A[i][j])) ? fabs(A[i][j]) : max[i];
        }
        if (fabs(max[i]) < 1e-15)
        {
            fprintf(stderr, "GaussEliPPP: A is singular\n");
            return NULL;
        }
    }

    // do the gauss elimination with partial pivoting proportionally
    for (int k = 0; k < N - 1; k++)
    {
        int r = k;
        // choose the pivot element
        for (int i = k; i < N; i++)
        {
            r = (fabs(A[r][k] / max[r]) < fabs(A[i][k] / max[i])) ? i : r;
        }
        if (fabs(A[r][k]) < 1e-15)
        {
            printf("A is singular\n");
            return NULL;
        }
        if (r != k)
        {
            double q = max[k]; max[k] = max[r]; max[r] = q;
            double *t = A[k]; A[k] = A[r]; A[r] = t;
        }

        for (int i = k + 1; i < N; i++)
        {
            A[i][k] = A[i][k] / A[k][k];
            for (int j = k + 1; j < N + 1; j++)
            {
                A[i][j] -= A[i][k] * A[k][j];
            }
        }
    }

    if (fabs(A[N - 1][N - 1]) < 1e-15)
    {
        printf("A is singular\n");
        return NULL;
    }
    
    for (int k = N - 1; k >= 0; k--) 
    {
        double t = 0;
        for (int j = k + 1; j < N; j++) 
        {
            t += A[k][j] * x[j];
        }
        x[k] = (A[k][N] - t) / A[k][k];
    }
    free(max);
    return x;
}

