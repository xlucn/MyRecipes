#include "NR.h"
#include "NRprivate.h"
/**
 * @brief Classic Runge-Kutta Method to solve a System of ODEs.
 * @param m the number of ODEs,
 * @param N the numebr of steps,
 * @param (a, b) the interval,
 * @param y0 the initial values,
 * @param f point to an array of functions.
 * @return a 2D array of values of all functions in all steps.
 */
double **SODERungeKutta(double (**f)(double, double*), double a, double b, double *y0, int m, int N)
{
    double h = (b - a) / N;
    double t = a;
    double *w = newArray1d(m);
    double **y = newArray2d(N + 1, m);
    double **k = newArray2d(4, m);

    for(int i = 0; i < m; i++)
    {
        w[i] = y0[i];
        y[0][i] = w[i];
    }
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < m; j++)
        {
            k[0][j] = h * f[j](t, w);
            w[j] = y[i][j] + k[0][j] / 2;
        }
        for(int j = 0; j < m; j++)
        {
            k[1][j] = h * f[j](t + h / 2, w);
            w[j] = y[i][j] + k[1][j] / 2;
        }
        for(int j = 0; j < m; j++)
        {
            k[2][j] = h * f[j](t + h / 2, w);
            w[j] = y[i][j] + k[2][j];
        }
        for(int j = 0; j < m; j++)
        {
            k[3][j] = h * f[j](t + h, w);
            y[i + 1][j] = y[i][j] + (k[0][j] + k[1][j] * 2 + k[2][j] * 2 + k[3][j]) / 6;
        }
        t += h;
    }
    delArray1d(w);
    delArray2d(4, N);
    return y;
}


