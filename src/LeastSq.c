#include <math.h>
#include <stdlib.h>
#include <NumericalRecipes.h>
#include <LibFunction.h>

/*
solve the leastsquare solution for a system of linear equations.
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
