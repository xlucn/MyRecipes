#include <math.h>
#include <stdio.h>
#include "NR.h"

/**
 * @brief Gauss elimination with partial pivoting
 * @param N The numbers of the variables
 * @param a The coefficient matrix in a 1-dimension array
 * @param b The constant vector
 * @return An array of numbers, NULL if the equation has no solution
 */
double *GaussEliPP(int N, double *a, double *b)
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
    
    int max;
    double *x = (double*)malloc_s(N * sizeof(double));

    for (int k = 0; k < N - 1; k++)
    {
        /* find the pivot */
        max = k;
        for (int i = k; i < N; i++)
        {
            if (fabs(A[i][k]) > fabs(A[max][k]))
            {
                max = i;
            }
        }
        /* singular matrix */
        if (fabs(A[max][k]) < 1e-15)
        {
            fprintf(stderr, "GaussEliPP: A is singular\n");
            return NULL;
        }
        /* swap */
        if (max != k) 
        {
            double *temp = A[k];
            A[k] = A[max];
            A[max] = temp;
        }
        /* elimination */
        for (int i = k + 1; i < N; i++) 
        {
            A[i][k] = A[i][k] / A[k][k];
            for (int j = k + 1; j < N + 1; j++) 
            {
                A[i][j] -= A[i][k] * A[k][j];
            }
        }
    }
    /* singular matrix */
    if (fabs(A[N - 1][N - 1]) < 1e-15) 
    {
        fprintf(stderr, "GaussEliPP: A is singular\n");
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
    return x;
}
