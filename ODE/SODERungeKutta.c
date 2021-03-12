/**
 * @file SODERungeKutta.c
 * @brief SODERungeKutta.c
 */

#include <stdlib.h>
#include "NR.h"

/**
 * @brief Classic Runge-Kutta Method to solve a System of ODEs.
 * @param f a pointer to an array of functions.
 * @param a lower limit of interval
 * @param b upper limit of interval
 * @param y0 the initial values,
 * @param m the number of ODEs,
 * @param N the numebr of steps,
 * @return an ::SODEsol variable
 */
SODEsol SODERungeKutta(double *(*f)(double, double*), double a, double b, double *y0, int m, int N)
{
    double h = (b - a) / N;
    double *t = malloc((N + 1) * sizeof(double));
    double *w = malloc(m * sizeof(double));
    double **y = newArray2d(N + 1, m);
    double **k = malloc(4 * sizeof(double*));

    t[0] = a;
    for(int i = 0; i < m; i++)
    {
        w[i] = y0[i];
        y[0][i] = w[i];
    }
    for(int i = 0; i < N; i++)
    {
        k[0] = f(t[i], w);

        for(int j = 0; j < m; j++)
        {
            w[j] = y[i][j] + k[0][j] * h / 2;
        }
        k[1] = f(t[i] + h / 2, w);

        for(int j = 0; j < m; j++)
        {
            w[j] = y[i][j] + k[1][j] * h / 2;
        }
        k[2] = f(t[i] + h / 2, w);

        for(int j = 0; j < m; j++)
        {
            w[j] = y[i][j] + k[2][j] * h;
        }
        k[3] = f(t[i] + h, w);

        for(int j = 0; j < m; j++)
        {
            y[i + 1][j] = y[i][j] + (k[0][j] + k[1][j] * 2 + k[2][j] * 2 + k[3][j]) * h / 6;
        }

        t[i + 1] = t[i] + h;
        free(k[0]);
        free(k[1]);
        free(k[2]);
        free(k[3]);
    }
    free(w);
    free(k);
    return newSODEsol(N + 1, t, y);
}


