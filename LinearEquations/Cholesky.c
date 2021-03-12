/**
 * @file Cholesky.c
 * @brief Cholesky.c
 */

#include <math.h>
#include "NR.h"

/**
 * @brief Decomposite a matrix with Cholesky method
 * @param N the numbers of the variables
 * @param a The symmetric matrix in a 1-dimension array
 */
void Cholesky(int N, double *a)
{
    for (int j = 0; j < N; j++)
    {
        /* calculate ajj */
        for (int k = 0; k < j; k++)
        {
            a[j * N + j] -= a[j * N + k] * a[j * N + k];
        }
        a[j * N + j] = sqrtf(a[j * N + j]);
        /* calculate aij, i >= j */
        for (int i = j + 1; i < N; i++)
        {
            for (int k = 0; k < j; k++)
            {
                a[i * N + j] -= a[i * N + k] * a[j * N + k];
            }
            a[i * N + j] /= a[j * N + j];
        }
    }
}
