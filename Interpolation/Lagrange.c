/**
 * @file Lagrange.c
 * @brief Lagrange.c
 */

#include <stdlib.h>
#include "NR.h"

/**
 * @brief Lagrange polynomial
 * @param a interpolation points
 * @param x variable value
 * @param N the number of intervals
 * @returns Lagrange polynomial value at x
 */
double *LagrangePoly(double *a, double x, int N)
{
    double *l = malloc(N * sizeof(double));
    for(int i = 0; i < N; i++)
    {
        l[i] = 1;
        for (int j = 0; j < N; j++)
        {
            if (j != i)
            {
                l[i] = l[i] * (x - a[j]) / (a[i] - a[j]);
            }
        }
    }
    return l;
}

