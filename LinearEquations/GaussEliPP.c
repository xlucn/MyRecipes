
#include <math.h>
#include <stdio.h>
#include "NR.h"

/**
 * @brief Gauss elimination with partial pivoting
 * @param N the rank of the equation
 * @param a the augmented matrix
 */
double *GaussEliPP(int N, double **a)
{
    int max;
    double *x = (double*)malloc_s(N * sizeof(double));

    for (int k = 0; k < N - 1; k++)
    {
        max = k;
        for (int i = k; i < N; i++)
        {
            if (fabs(a[i][k]) > fabs(a[max][k]))
            {
                max = i;
            }
        }
        if (a[max][k] == 0)
        {
            fprintf(stderr, "A is singular\n");
            return NULL;
        }
        if (max != k) 
        {
            double *temp = a[k];
            a[k] = a[max];
            a[max] = temp;
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
        fprintf(stderr, "A is singular\n");
        return NULL;
    }
    else
    {
        x[N - 1] = a[N - 1][N] / a[N - 1][N - 1];
    }
    for (int k = N - 2; k >= 0; k--) {
        double t = 0;
        for (int j = k + 1; j < N; j++) {
            t += a[k][j] * x[j];
        }
        x[k] = (a[k][N] - t) / a[k][k];
    }
    return x;
}
