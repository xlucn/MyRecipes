#include <math.h>
#include <stdio.h>
#include "NR.h"


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
        if (fabs(a[max[k]][k]) < FLOAT_ZERO_LIM)
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
