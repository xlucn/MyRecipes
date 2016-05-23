#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "NR.h"
#include "Test.h"

/**
 * testing orthogonal decomposition
 */

/**
 * @brief 
 * @returns 
 * 
 * 
 */
int testLeastSq()
{
    int m = 4;
    int n = 3;
    double temp[4][3] = {{1, -2, 1},{0, 1, -1},{2, -4, 3},{4, -7, 4}};
    double **A = (double**)malloc_s(m * sizeof(double));
    for (int i = 0; i < m; i++)
    {
        A[i] = temp[i];
    }
    double b[4] = {-4, 3, 1, -6};

    double *res = LeastSquare(m, n, A, b);
    for (int i = 0; i < n; i++)
    {
        printf("%lf\t", res[i]);
    }
    printf("\n");
    return PASSED;
}

/**
 * @brief 
 * @returns 
 * 
 * 
 */
int testQR()
{
    int m = 4;
    int n = 3;
    double eps = 1e-6;
    double temp[4][3] = {{1,-1,1},{0,-2,-1},{1,1,1},{0,2,1}};
    double **A = (double**)malloc_s(m * sizeof(double*));
    for(int i = 0; i < m; i++)
    {
        A[i] = temp[i];
    }

    double ***QR = ImprovedGramSchmidtQR(m, n, A);
    for(int i = 0; i < m; i++)
    {
        for(int j = 0; j < n; j++)
        {
            printf("%lf\t",QR[0][i][j]);
        }
        printf("\n");
    }
    printf("\n");
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            printf("%lf\t",QR[1][i][j]);
        }
        printf("\n");
    }
    printf("\n");
    for(int i = 0; i < m; i ++)
    {
        for(int j = 0; j < n; j++)
        {
            for(int k = 0; k < n; k++)
            {
                A[i][j] -= QR[0][i][k] * QR[1][k][j];
            }
            if(fabs(A[i][j]) > eps)
            {
                return FAILED;
            }
        }
    }
    return PASSED;
}


