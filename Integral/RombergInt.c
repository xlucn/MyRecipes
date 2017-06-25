/** @file RombergInt.c */
#include <math.h>
#include "NR.h"

/**
 * @brief Romberg integration
 * @param f the integration function
 * @param a the start point
 * @param b the end point of interval,
 * @param N max steps,
 * @param eps the precision,
 * @return the integral
 */
double RombergInt(double (*f)(double), double a, double b, int N, double eps)
{
    double h = b - a;
    double *T1 = newArray1d(1), *T2;
    T1[0] = (f(b) + f(a)) * h / 2;

    for (int i = 1; i < N; i++, h /= 2)
    {
        T2 = newArray1d(i + 1);
        //T[i][0]
        T2[0] = T1[0];
        for (int k = 0; k < pow(2, i - 1) - 0.5; k++)
        {
            T2[0] += h * f(a + (k + 0.5) * h);
        }
        T2[0] /= 2;

        //T[i][j]
        for (int j = 1; j < i + 1; j++)
        {
            T2[j] = (pow(4, j) * T2[j - 1] - T1[j - 1]) / (pow(4, j) - 1);
        }

        if (fabs(T2[i] - T2[i - 1]) <= eps)
        {
            double result = T2[i];
            delArray1d(T1);
            delArray1d(T2);
            return result;
        }

        delArray1d(T1);
        T1 = T2;
    }
    double result = T1[N - 1];
    delArray1d(T1);
    return result;
}
