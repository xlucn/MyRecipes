/** @file GramSchmidtQR.c */
#include <math.h>
#include "NR.h"
#include "NRprivate.h"

/**
 * @brief Gram-Schmidt method of POD(proper orthogonal decomposition).
 * @param m length of first axis of A
 * @param n length of second axis of A
 * @param A matrix of shape m * n (m > n) and rank(A) = n
 * @returns a three demension array which contains two 2-d arrays [Q,R].
 */
double ***GramSchmidtQR(int m, int n, double **A)
{
    double **Q = newArray2d(m, n);
    double **R = newArray2d(n, n);
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
