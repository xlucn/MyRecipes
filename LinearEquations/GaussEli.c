#include <math.h>
#include "NR.h"

/**
 * @brief Gaussian elimination method to solve linear equation in the form of Ax=b.
 * @param N The numbers of the variables
 * @param a The coefficient matrix in a 1-dimension array
 * @param b The constant vector
 * @returns An array of numbers, NULL if the equation has no solution
 */
double *GaussEli(int N, double *a, double *b)
{
    double **A = AugmentedMatrix(a, b, N, N, 1);
    double *x = newArray1d(N);
    double l;

    for(int k = 0; k < N - 1; k++)
    {
        /* if a_kk is zero, pick another nonzero from this row */
        if(fabs(A[k][k]) < FLOAT_ZERO_LIM)
        {
            for(int i = k + 1; i < N; i++)
            {
                if(fabs(A[i][k]) >= FLOAT_ZERO_LIM)
                {
                    double *temp = A[i];
                    A[i] = A[k];
                    A[k] = temp;
                    break;
                }
            }
            if(fabs(A[k][k]) < FLOAT_ZERO_LIM)
            {
                return NULL;
            }
        }
        /* elimination */
        for(int i = k + 1; i < N; i++)
        {
            l = A[i][k] / A[k][k];
            A[i][k] = 0;
            for(int j = k + 1; j < N + 1; j++)
            {
                A[i][j] -= l * A[k][j];
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

    delArray2d(A, N);
    return x;
}
