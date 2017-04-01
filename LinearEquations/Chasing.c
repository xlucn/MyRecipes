#include "NR.h"


/**
 * @brief Chasing method for solving tridiagonal equations.
 * 
 * - d0  c0   0  ...  0  -     - b1 -
 * | a1  d1  c1  ...  0  |     | b2 |
 * |  0  a2  d2  ... ... | x = | .. |
 * | ... ... ... ... cN-1|     | .. |
 * -  0   0  ... aN  dN  -     - bn -
 * 
 * @param N the order of the matrix.
 * @param d the main diagonal,
 * @param c the line above d,
 * @param a the line below d,
 * @param b the constant vector.
 * @returns the solution of the equation *x or NULL if the eqation has no solution.
 */
double *Chasing(int N, double *d, double *c, double *a, double *b)
{
    double *p = (double *)malloc_s(N * sizeof(double));
    double *q = (double *)malloc_s(N * sizeof(double));
    double *x = (double *)malloc_s(N * sizeof(double));
    double *y = (double *)malloc_s(N * sizeof(double));

    if(d[0] == 0)
    {
        return NULL;
    }
    p[0] = d[0];
    q[0] = c[0] / d[0];
    for(int i = 1; i < N - 1; i++)
    {
        p[i] = d[i] - a[i] * q[i - 1];
        if(p[i] == 0)
        {
            return NULL;
        }
        q[i] = c[i] / p[i];
    }
    p[N - 1] = d[N - 1] - a[N - 1] * q[N - 2];
    if(p[N - 1] == 0)
    {
        return NULL;
    }
    y[0] = b[0] / p[0];
    for(int i = 1; i < N; i++)
    {
        y[i] = (b[i] - a[i] * y[i - 1]) / p[i];
    }
    x[N - 1] = y[N - 1];
    for(int i = N - 2; i >= 0; i--)
    {
        x[i] = y[i] - q[i] * x[i + 1];
    }
    free(p);
    free(q);
    free(y);
    return x;
}
