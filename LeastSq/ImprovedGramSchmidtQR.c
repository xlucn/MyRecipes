/** @file ImprovedGramSchmidtQR.c */
#include <math.h>
#include "NRprivate.h"
#include "NR.h"

/**
 * @brief Improved version of Gram-Schmidt method.
 * @param m length of first axis of A
 * @param n length of second axis of A
 * @param A matrix of shape m * n (m > n) and rank(A) = n
 * @returns a three demension array which contains two 2-d arrays [Q,R].
 */
double ***ImprovedGramSchmidtQR(int m, int n, double **A)
{
    double temp;
    double ***QR = (double***)malloc_s(2 * sizeof(double**));
    double **Q = newArray2d(m, n);
    double **R = newArray2d(n, n);
    QR[0] = Q;
    QR[1] = R;
    for(int i = 0; i < m; i ++)
    {
        for(int j = 0; j < n; j++)
        {
            Q[i][j] = A[i][j];
        }
    }

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
