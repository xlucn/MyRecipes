#include <math.h>
#include <stdlib.h>
#include "NR.h"

/**
 * @brief solve the leastsquare solution for a system of linear equations.
 */
double *LeastSquare(int m, int n, double **A, double *b)
{
    if (m<n)
    {
        return NULL;
    }

    double *res = (double*)malloc_s(n * sizeof(double));
    double *Ab = (double*)malloc_s(n * sizeof(double));
    double **AA = (double**)malloc_s(n * sizeof(double*));
    for (int i = 0; i < n; i++)
    {
        AA[i] = (double*)malloc_s(n * sizeof(double));
    }

    //transpose of matrix A times A
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            AA[i][j] = 0;
            for (int k = 0; k < m; k++)
            {
                AA[i][j] += A[k][i] * A[k][j];
            }
        }
    }

    //transpose of matrix A times b
    for (int i = 0; i < n; i++)
    {
        Ab[i] = 0;
        for (int j = 0; j < m; j++)
        {
            Ab[i] += A[j][i] * b[j];
        }
    }

    //solve the equation AA x = Ab
    res = GaussEli(n, AA, Ab);

    //free the memory space allocated before
    for (int i = 0; i < n; i++)
    {
        free(AA[i]);
    }
    free(AA);
    free(Ab);
    return res;
}
