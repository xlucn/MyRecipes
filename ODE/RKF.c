/** @file RKF.c */
#include <math.h>
#include <stdio.h>
#include "NR.h"
#include "NRprivate.h"
#include "constants.h"
/**
 * @brief RKF method(Ronge Kutta Fehlberg method)
 * @param f right function of ODE
 * @param a lower limit of interval
 * @param b upper limit of interval
 * @param y0 initial value of f
 * @param TOL tolerance
 * @param hmax max step size
 * @param hmin min step size
 * @param A the paramaters in A matrix of RKF method
 * @param B the paramaters in B vector of RKF method
 * @param Bstar the paramaters in Bstar vector of RKF method
 * @param C the paramaters in C vector of RKF method
 * @param n the size of A matrix, or the length of B, Bstar and C
 * @returns integral
 */
static ODEsol RKFmn(double (*f)(double,double), double a, double b, double y0, double TOL, double hmax, double hmin,
    const double* A, const double* B, const double* Bstar, const double* C, int n)
{
    int step = 0; //the total steps
    int length = 1024;
    double t = a;
    double y = y0;
    double h = hmax;
    double *k = newArray1d(n);
    double R = 0;
    double delta;
    step++;
    double *ts = newArray1d(length);
    double *ys = newArray1d(length);
    ts[0] = t;
    ys[0] = y;

    while (t < b)
    {
        for(int i = 0; i < n; i++)
        {
            y = ys[step - 1];
            for(int j = 0; j < i; j++)
            {
                y += A[i * n + j] * k[j];
            }
            k[i] = h * f(t + C[i] * h, y);
            R += (B[i] - Bstar[i]) * k[i] / h;
        }
        y = ys[step - 1];
        R = fabs(R);
        delta = pow((TOL / R / 2.0), 0.25);

        if (R <= TOL)
        {
            t += h;
            for(int i = 0; i < n; i++)
            {
                y += B[i] * k[i];
            }
            step++;
            if(step == length)
            {
                length += 1024;
                ts = Array1dResize(ts, length);
                ys = Array1dResize(ys, length);
            }
            ts[step - 1] = t;
            ys[step - 1] = y;
        }

        h = delta < 0.1 ? 0.1 * h : (delta > 4 ? 4 * h: delta * h);

        if (h >= hmax)
        {
            h = hmax;
        }
        if (h < hmin)
        {
            printf("h is smaller than hmin\n");
            return NULL;
        }
    }
    ODEsol sol = newODEsol(step, ts, ys);
    return sol;
}

/**
 * @brief RKF method(Ronge Kutta Fehlberg method) of orders 8 and 7
 * @param f right function of ODE
 * @param a lower limit of interval
 * @param b upper limit of interval
 * @param y0 initial value of f
 * @param TOL tolerance
 * @param hmax max step size
 * @param hmin min step size
 * @returns the integral
 */
ODEsol RKF78(double(*f)(double,double), double a, double b, double y0, double TOL, double hmax, double hmin)
{
    return RKFmn(f, a, b, y0, TOL, hmax, hmin, A78, B78, Bstar78, C78, 13);
}

/**
 * @brief RKF method(Ronge Kutta Fehlberg method) with orders 5 and 4
 * @param f right function of ODE
 * @param a lower limit of interval
 * @param b upper limit of interval
 * @param y0 initial value of f
 * @param TOL tolerance
 * @param hmax max step size
 * @param hmin min step size
 * @returns the integral
 */
ODEsol RKF45(double(*f)(double,double), double a, double b, double y0, double TOL, double hmax, double hmin)
{
    return RKFmn(f, a, b, y0, TOL, hmax, hmin, A45, B45, Bstar45, C45, 6);
}

