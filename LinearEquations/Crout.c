#include <math.h>
#include "NR.h"
#include "NRprivate.h"
#include "constants.h"
/**
 * @brief Solve a system of linear equations with Crout method
 * @param N the numbers of the variables
 * @param a The coefficient matrix in a 1-dimension array
 * @param b The constant vector
 * @returns An array of numbers, NULL if a is singular
 */
double *Crout(int N, double *a, double *b)
{
    double **LU = AugmentedMatrix(a, b, N, N, 1);
    double *x = newArray1d(N);
    
    for(int k = 0; k < N; k++)
    {
        /* the kth col of L */
        for(int i = k; i < N; i++)
        {
            for(int r = 0; r < k; r++)
            {
                LU[i][k] -= LU[i][r] * LU[r][k];
            }
        }
        
        /* choose the pivot */
        int pivot = k;
        for(int i = k; i < N; i++)
        {
            if(fabs(LU[i][k]) > fabs(LU[pivot][k]))
            {
                pivot = i;
            }
        }
        if(fabs(LU[pivot][k]) < FLOAT_ZERO_LIM)
        {
            return NULL;
        }
        if (pivot != k)
        {
            double *temp = LU[k];
            LU[k] = LU[pivot];
            LU[pivot] = temp;
        }
        
        /* the kth row of U and solve Ly=b for y at the same time*/
        for(int j = k + 1; j <= N; j++)
        {
            for(int r = 0; r < k; r ++)
            {
                LU[k][j] -= LU[k][r] * LU[r][j];
            }
            LU[k][j] /= LU[k][k];
        }
    }
    
    /* solve Ux=y for x */
    for(int k = N - 1; k >= 0; k--)
    {
        x[k] = LU[k][N];
        for(int r = k + 1; r < N; r++)
        {
            x[k] -= LU[k][r] * x[r];
        }
    }
    
    delArray2d(LU, N);
    return x;
}
