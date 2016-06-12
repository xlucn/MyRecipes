#include <math.h>
#include <stdio.h>
#include "NR.h"

/**
 * @brief Gauss elimination with partial pivoting proportionally
 * @param N the rank of the equation
 * @param a the augmented matrix
 * @return an array of numbers, NULL if the equation has no solution
 */
double *GaussEliPPP(int N, double **a)
{
    double *max = (double*)malloc_s(N * sizeof(double));
    double *x = (double*)malloc_s(N * sizeof(double));

    /* get the largest number of each line */
    for (int i = 0; i < N; i++)
    {
        max[i] = 0;
        for (int j = 0; j < N; j++)
        {
            max[i] = (fabs(max[i]) < fabs(a[i][j])) ? fabs(a[i][j]) : max[i];
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
            r = (fabs(a[r][k] / max[r]) < fabs(a[i][k] / max[i])) ? i : r;
        }
        if (fabs(a[r][k]) < 1e-15)
        {
            printf("A is singular\n");
            return NULL;
        }
        if (r != k)
        {
            double q = max[k]; max[k] = max[r]; max[r] = q;
            double *t = a[k]; a[k] = a[r]; a[r] = t;
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

    if (fabs(a[N - 1][N - 1]) < 1e-15)
    {
        printf("A is singular\n");
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
    free(max);
    return x;
}

