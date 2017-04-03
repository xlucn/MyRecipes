#include <math.h>
#include <stdio.h>
#include "NR.h"
#include "constants.h"
/**
 * @brief 
 * @param double 
 * @param double 
 * @param a 
 * @param b 
 * @param y0 
 * @param TOL 
 * @param hmax 
 * @param hmin 
 * @param A 
 * @param B 
 * @param Bstar 
 * @param C 
 * @param n 
 * @returns 
 */
static double *RKFmn(double(*f)(double,double), double a, double b, double y0, double TOL, double hmax, double hmin,
    const double* A, const double* B, const double* Bstar, const double* C, int n)
{
    int step = 0; //the total steps
    double t = a;
    double y = y0;
    double h = hmax;
    double *k = newArray1d(n);
    double R = 0;
    double delta;
    double *result;
    step++;
    result = newArray1d(step * 3 + 1);
    result[1] = t;
    result[2] = h;
    result[3] = y;

    while (t < b)
    {
        for(int i = 0; i < n; i++)
        {
            y = result[3 * step];
            for(int j = 0; j < i; j++)
            {
                y += A[i * n + j] * k[j];
            }
            k[i] = h * f(t + C[i] * h, y);
            R += (B[i] - Bstar[i]) * k[i] / h;
        }
        y = result[3 * step];
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
            result = (double*)realloc_s(result, (step * 3 + 1) * sizeof(double));
            result[step * 3 - 2] = t;
            result[step * 3 - 1] = h;
            result[step * 3] = y;
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
    result[0] = step;
    return result;
}
/**
 * @brief Runge-Kutta-Fehlberg Method
 */
double *RKF78(double(*f)(double,double), double a, double b, double y0, double TOL, double hmax, double hmin)
{
    return RKFmn(f, a, b, y0, TOL, hmax, hmin, A78, B78, Bstar78, C78, 13);
}

double *RKF45(double(*f)(double,double), double a, double b, double y0, double TOL, double hmax, double hmin)
{
    return RKFmn(f, a, b, y0, TOL, hmax, hmin, A45, B45, Bstar45, C45, 6);
}

