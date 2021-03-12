/**
 * @file GaussJordanEli.c
 * @brief GaussJordanEli.c
 */

#include <stdlib.h>
#include <math.h>
#include "NR.h"
#include "constants.h"

/**
 * @brief Gauss Jordan elimination method to solve a system of linear equations
 * @param N the rank of the equation
 * @param A the coefficient matrix in 1d form
 * @param b the constant vector
 * @returns An array of numbers, NULL if the equation has no solution
 */
double *GaussJordanEli(int N, double *A, double *b)
{
    double **Ab = AugmentedMatrix(A, b, N, N, 1);
    double* x = malloc(N * sizeof(double));

    for (int k = 0; k < N; k++)
    {
        int pivot = k;
        for (int i = k; i < N; i++) if(fabs(Ab[pivot][k]) < fabs(Ab[i][k]))
        {
            pivot = i;
        }
        double *temp = Ab[k];
        Ab[k] = Ab[pivot];
        Ab[pivot] = temp;
        if (fabs(Ab[k][k]) < FLOAT_ZERO_LIM)
        {
            return NULL;
        }

        for(int j = k + 1; j <= N; j++)
        {
            Ab[k][j] /= Ab[k][k];
            for(int i = 0; i < N; i++) if(i != k)
            {
                Ab[i][j] -= Ab[i][k] * Ab[k][j];
            }
        }
    }

    for (int k = 0; k < N; k++)
    {
        x[k] = Ab[k][N];
    }

    delArray2d(Ab, N);
    return x;
}
