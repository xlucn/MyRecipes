#include <math.h>
#include <stdlib.h>
#include "NR.h"

/**
 * @brief Gram-Schmidt method of POD(proper orthogonal decomposition).
 * Rank of The m*n(m>n) matrix need to be n.
 * return a three demension array which contains two 2-d arrays [Q,R].
 */
double ***GramSchmidtQR(int m, int n, double **A)
{
    double **Q = (double**)malloc_s(m * sizeof(double*));
    for(int i = 0; i < m; i ++)
    {
        Q[i] = (double*)malloc_s(n * sizeof(double));
    }
    double **R = (double**)malloc_s(n * sizeof(double*));
    for(int i = 0; i < n; i ++)
    {
        R[i] = (double*)malloc_s(n * sizeof(double));
    }
    double ***QR = (double***)malloc_s(2 * sizeof(double**));
    QR[0] = Q;
    QR[1] = R;

    for(int j = 0; j < n; j++)
    {
        R[j][j] = 0;
        for(int k = 0; k < m; k++)
        {
            Q[k][j] = 0;
        }
        for(int i = 0; i < j; i++)
        {
            R[i][j] = 0;
            for(int k = 0; k < m; k++)
            {
                R[i][j] += Q[k][i] * A[k][j];
            }
            for(int k = 0; k < m; k++)
            {
                Q[k][j] += R[i][j] * Q[k][i];
            }
        }
        for(int k = 0; k < m; k++)
        {
            Q[k][j] = A[k][j] - Q[k][j];
            R[j][j] += Q[k][j] * Q[k][j];
        }
        R[j][j] = sqrt(R[j][j]);
        for(int k = 0; k < m; k++)
        {
            Q[k][j] = Q[k][j] / R[j][j];
        }
    }
    return QR;
}
