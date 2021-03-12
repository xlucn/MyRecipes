/**
 * @file GaussEliPP.c
 * @brief GaussEliPP.c
 */

#include <stdlib.h>
#include <math.h>
#include "NR.h"
#include "constants.h"

/**
 * @brief Gauss elimination with partial pivoting
 * @param N The numbers of the variables
 * @param a The coefficient matrix in a 1-dimension array
 * @param b The constant vector
 * @return An array of numbers, NULL if the equation has no solution
 */
double *GaussEliPP(int N, double *a, double *b)
{
    double **A = AugmentedMatrix(a, b, N, N, 1);
    double *x = malloc(N * sizeof(double));

    for (int k = 0; k < N; k++)
    {
        /* find the pivot */
        int max = k;
        for (int i = k; i < N; i++)
        {
            if (fabs(A[i][k]) > fabs(A[max][k]))
            {
                max = i;
            }
        }
        /* singular matrix */
        if (fabs(A[max][k]) < FLOAT_ZERO_LIM)
        {
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

    for (int k = N - 1; k >= 0; k--)
    {
        double t = 0;
        for (int j = k + 1; j < N; j++)
        {
            t += A[k][j] * x[j];
        }
        x[k] = (A[k][N] - t) / A[k][k];
    }

    delArray2d(A, N);
    return x;
}
