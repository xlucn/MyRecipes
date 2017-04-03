#include <math.h>
#include "NR.h"

/**
 * @brief Gauss Jordan elimination method to solve a system of linear equations
 * @param N the rank of the equation
 * @param a the augmented matrix
 */
double *GaussJordanEli(int N, double *A, double *b)
{
    double **Ab = AugmentedMatrix(A, b, N, N, 1);
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
                max[k] = (fabs(Ab[max[k]][k]) < fabs(Ab[i][k])) ? i : max[k];
            }
        }
        if (fabs(Ab[max[k]][k]) < FLOAT_ZERO_LIM)
        {
            return NULL;
        }

        for (int i = 0; i < N; i++)
        {
            if (i != max[k])
            {
                Ab[i][k] = Ab[i][k] / Ab[max[k]][k];
                for (int j = k + 1; j < N + 1; j++)
                {
                    Ab[i][j] -= Ab[i][k] * Ab[max[k]][j];
                }
            }

        }
    }

    for (int k = 0; k < N; k++)
    {
        x[k] = Ab[max[k]][N] / Ab[max[k]][k];
    }
    
    free(max);
    delArray2d(Ab, N);
    return x;
}
