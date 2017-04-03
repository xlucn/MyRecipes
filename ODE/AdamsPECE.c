#include "NR.h"
/**
 * @brief One-step Adams correlation PECE method. Use classic Runge-Kutta method for the initial value.
 * @param f right function of ODE
 * @param a lower limit of interval
 * @param b upper limit of interval
 * @param dy0 initial value of derivative of f
 * @param y0 initial value of f
 * @param N number of subintervals
 * @return array of ys
 */
double *AdamsPECE(double(*f)(double, double), double a, double b, double dy0, double y0, int N)
{
    double *y = newArray1d(N + 1);
    double *dy = newArray1d(N + 1);
    double h = (b - a) / N;
    y[0] = y0;
    y[1] = ClassicRungeKutta(f, a, b, y0, N)[1];
    dy[0] = dy0;
    dy[1] = f(a + h, y[1]);

    for(int i = 1; i < N; i ++)
    {
        y[i + 1] = y[i] + h / 2 * (3 * dy[i] - dy[i - 1]);
        dy[i + 1] = f(a + (i + 1) * h, y[i + 1]);
        y[i + 1] = y[i] + h / 2 * (dy[i] + dy[i + 1]);
        dy[i + 1] = f(a + (i + 1) * h, y[i + 1]);
    }
    delArray1d(dy);
    return y;
}

