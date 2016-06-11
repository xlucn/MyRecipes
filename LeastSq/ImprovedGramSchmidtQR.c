#include <math.h>
#include <stdlib.h>
#include "NR.h"

/**
 * @brief Improved version of Gram-Schmidt method.
 */
double ***ImprovedGramSchmidtQR(int m, int n, double **A)
{
    double temp;
    double ***QR = (double***)malloc_s(2 * sizeof(double**));
    double **Q = (double**)malloc_s(m * sizeof(double*));
    for(int i = 0; i < m; i ++)
    {
        Q[i] = (double*)malloc_s(n * sizeof(double));
        for(int j = 0; j < n; j++)
        {
            Q[i][j] = A[i][j];
        }
    }
    double **R = (double**)malloc_s(n * sizeof(double*));
    for(int i = 0; i < n; i ++)
    {
        R[i] = (double*)malloc_s(n * sizeof(double));
    }
    QR[0] = Q;
    QR[1] = R;

    for(int k = 0; k < n; k++)
    {
        temp = 0;
        for(int i = 0; i < m; i++)
        {
            temp += Q[i][k] * Q[i][k];
        }
        R[k][k] = sqrt(temp);
        for(int i = 0; i < m; i++)
        {
            Q[i][k] = Q[i][k] / R[k][k];
        }
        for(int j = k + 1; j < n; j++)
        {
            R[k][j] = 0;
            for(int l = 0; l < m; l++)
            {
                R[k][j] += Q[l][k] * Q[l][j];
            }
            for(int l = 0; l < m; l++)
            {
                Q[l][j] -= R[k][j] * Q[l][k];
            }
        }
    }
    return QR;
}
