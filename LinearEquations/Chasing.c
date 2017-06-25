/** @file Chasing.c */
#include <math.h>
#include "constants.h"
#include "NR.h"
#include "NRprivate.h"


/**
 * @brief Chasing method for solving tridiagonal equations.
 * @param N the order of the matrix.
 * @param d the main diagonal,
 * @param c the line above d,
 * @param a the line below d,
 * @param b the constant vector.
 * @returns the solution of the equation *x or NULL if the eqation has no solution.
 */
double *Chasing(int N, double *d, double *c, double *a, double *b)
{
    double *p = newArray1d(N);
    double *q = newArray1d(N);
    double *x = newArray1d(N);
    double *y = newArray1d(N);

    if(fabs(d[0]) < FLOAT_ZERO_LIM){ return NULL; }
    p[0] = d[0];
    q[0] = c[0] / d[0];
    for(int i = 1; i < N - 1; i++)
    {
        p[i] = d[i] - a[i] * q[i - 1];
        if(fabs(p[i]) < FLOAT_ZERO_LIM) { return NULL; }
        q[i] = c[i] / p[i];
    }
    p[N - 1] = d[N - 1] - a[N - 1] * q[N - 2];
    if(fabs(p[N - 1]) < FLOAT_ZERO_LIM){ return NULL; }
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
    delArray1d(p);
    delArray1d(q);
    delArray1d(y);
    return x;
}
