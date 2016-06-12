#include <math.h>
#include <stdio.h>
#include "NR.h"

/**
 * @brief Gauss elimination with partial pivoting
 * @param N the rank of the equation
 * @param a the augmented matrix
 * @return an array of numbers, NULL if the equation has no solution
 */
double *GaussEliPP(int N, double **a)
{
    int max;
    double *x = (double*)malloc_s(N * sizeof(double));

    for (int k = 0; k < N - 1; k++)
    {
        /* find the pivot */
        max = k;
        for (int i = k; i < N; i++)
        {
            if (fabs(a[i][k]) > fabs(a[max][k]))
            {
                max = i;
            }
        }
        /* singular matrix */
        if (fabs(a[max][k]) < 1e-15)
        {
            fprintf(stderr, "GaussEliPP: A is singular\n");
            return NULL;
        }
        /* swap */
        if (max != k) 
        {
            double *temp = a[k];
            a[k] = a[max];
            a[max] = temp;
        }
        /* elimination */
        for (int i = k + 1; i < N; i++) 
        {
            a[i][k] = a[i][k] / a[k][k];
            for (int j = k + 1; j < N + 1; j++) 
            {
                a[i][j] -= a[i][k] * a[k][j];
            }
        }
    }
    /* singular matrix */
    if (fabs(a[N - 1][N - 1]) < 1e-15) 
    {
        fprintf(stderr, "GaussEliPP: A is singular\n");
        return NULL;
    }
    
    for (int k = N - 1; k >= 0; k--) 
    {
        double t = 0;
        for (int j = k + 1; j < N; j++) 
        {
            t += a[k][j] * x[j];
        }
        x[k] = (a[k][N] - t) / a[k][k];
    }
    return x;
}
