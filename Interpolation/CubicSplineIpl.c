/** @file CubicSplineIpl.c */
#include <stdlib.h>
#include <math.h>
#include "NR.h"

/**
 * @brief Cubic Spline Interpolation.
 * @param N the number of data points
 * @param f a function
 * @param x variable value
 * @param a data points
 * @param m cubic spline parameters
 * @returns the value of interpolation poly at x
 *
 * the parameter *m is created by functions LagrangeCubicSplineIplPara(),
 * CompleteCubicSplineIplPara() and NatureCubicSplineIplPara().
 */
static double CubicSplineIpl(int N, double (*f)(double), double x, double *a, double *m)
{
    int i = 0;
    for (int j = 0; j < N; j++)
        if (x < a[j + 1] && x >= a[j])
        {
            i = j;
            break;
        }

    double h = a[i + 1] - a[i];

    return 1 / h * (m[i] / 6 * pow(a[i + 1] - x, 3) + m[i + 1] / 6 * pow(x - a[i], 3)) +
           f(a[i]) + DividedDiff(f, a + i, 2) * (x - a[i]) -
           h * h / 6 * ((m[i + 1] - m[i]) * (x - a[i]) / h + m[i]);
}

/**
 * @brief Generate the parameters for Lagrange Cubic Spline Interpolation.
 * @param N the nubmer of intervals
 * @param a the array on interpolation points
 * @param f the interpolation function
 * @return the parameters for Lagrange Cubic Spline Interpolation.
 */
static double *LagrangeCubicSplineIplPara(int N, double *a, double (*f)(double))
{
    double *h = malloc(N * sizeof(double));			/* step */
    double *d = malloc((N - 1) * sizeof(double));		/* constant vector */
    double *Matd = malloc((N - 1) * sizeof(double));	/* diagonal */
    double *Matc = malloc((N - 1) * sizeof(double));	/* the line above diagonal */
    double *Mata = malloc((N - 1) * sizeof(double));	/* the line below diagonal */

    for (int i = 0; i < N; i++)
    {
        h[i] = a[i + 1] - a[i];
    }

    /* d */
    d[0] = DividedDiff(f, a, 3) - (h[0] + h[1]) * DividedDiff(f, a, 4);
    for (int i = 1; i < N; i++)
    {
        d[i] = DividedDiff(f, a + i, 2) - DividedDiff(f, a + i - 1, 2);
    }
    d[N] = DividedDiff(f, a + N - 2, 3) + (h[N - 1] + h[N - 2]) * DividedDiff(f, a + N - 3, 4);
    for (int i = 0; i < N + 1; i++)
    {
        d[i] = d[i] * 6;
    }

    /* Matd */
    Matd[0] = 2;
    for (int i = 1; i < N; i++)
    {
        Matd[i] = (h[i - 1] + h[i]) * 2;
    }
    Matd[N] = 2;

    /* Mata */
    Mata[0] = 0;
    for (int i = 1; i < N; i++)
    {
        Mata[i] = h[i - 1];
    }
    Mata[N] = 1;

    /* Matc */
    Matc[0] = 1;
    for (int i = 1; i < N; i++)
    {
        Matc[i] = h[i];
    }
    Matc[N] = 0;

    /* solve the tridiagonal equations */
    double *result = Chasing(N + 1, Matd, Matc, Mata, d);
    free(Mata);
    free(Matc);
    free(Matd);
    free(d);
    free(h);
    return result;
}
/**
 * @brief Generate the parameters for complete Cubic Spline Interpolation.
 * @param N the nubmer of intervals
 * @param a the array on interpolation points
 * @param f the interpolation function
 * @param df_a derivative at a
 * @param df_b derivative at b
 * @return the parameters for complete Cubic Spline Interpolation.
 */
static double *CompleteCubicSplineIplPara(int N, double *a, double (*f)(double), double df_a, double df_b)
{
    double *h = malloc(N * sizeof(double));			/* step */
    double *d = malloc((N - 1) * sizeof(double));		/* constant vector */
    double *Matd = malloc((N - 1) * sizeof(double));	/* diagonal */
    double *Matc = malloc((N - 1) * sizeof(double));	/* the line above diagonal */
    double *Mata = malloc((N - 1) * sizeof(double));	/* the line below diagonal */

    for (int i = 0; i < N; i++)
    {
        h[i] = a[i + 1] - a[i];
    }

    Matd[0] = 2 * h[0];
    for (int i = 1; i < N; i++)
    {
        Matd[i] = 2 * (h[i - 1] + h[i]);
    }
    Matd[N] = 2 * h[N - 1];

    Mata[0] = 0;
    for (int i = 1; i < N + 1; i++)
    {
        Mata[i] = h[i - 1];
    }

    for (int i = 0; i < N; i++)
    {
        Matc[i] = h[i];
    }
    Matc[N] = 0;

    d[0] = DividedDiff(f, a, 2) - df_a;
    for (int i = 1; i < N; i++)
    {
        d[i] = DividedDiff(f, a + i, 2) - DividedDiff(f, a + i - 1, 2);
    }
    d[N] = df_b - DividedDiff(f, a + N - 1, 2);
    for (int i = 0; i < N + 1; i++)
    {
        d[i] = d[i] * 6;
    }

    /* solve the tridiagonal equations */
    double *result = Chasing(N + 1, Matd, Matc, Mata, d);
    free(Mata);
    free(Matc);
    free(Matd);
    free(d);
    free(h);
    return result;
}

/**
 * @brief Generate the parameters for nature Cubic Spline Interpolation.
 * @param N the nubmer of intervals
 * @param a the array on interpolation points
 * @param f the interpolation function
 * @param ddf_a second order derivative at a
 * @param ddf_b second order derivative at b
 * @return the parameters for nature Cubic Spline Interpolation.
 */
static double *NatureCubicSplineIplPara(int N, double *a, double (*f)(double), double ddf_a, double ddf_b)
{
    double *m = malloc((N + 1) * sizeof(double));		/* interpolation params */
    double *temp = malloc((N - 1) * sizeof(double));
    double *h = malloc(N * sizeof(double));			/* step */
    double *d = malloc((N - 1) * sizeof(double));		/* constant vector */
    double *Matd = malloc((N - 1) * sizeof(double));	/* diagonal */
    double *Matc = malloc((N - 1) * sizeof(double));	/* the line above diagonal */
    double *Mata = malloc((N - 1) * sizeof(double));	/* the line below diagonal */
    for (int i = 0; i < N; i++)
    {
        h[i] = a[i + 1] - a[i];
    }

    for (int i = 0; i < N - 1; i++)
    {
        Matd[i] = 2 * (h[i] + h[i + 1]);
        Mata[i] = h[i];
        Matc[i] = h[i + 1];
    }

    for (int i = 0; i < N - 1; i++)
    {
        d[i] = (DividedDiff(f, a + i + 1, 2) - DividedDiff(f, a + i, 2)) * 6;
    }
    d[0] = d[0] - h[0] * ddf_a;
    d[N - 2] = d[N - 2] - h[N - 1] * ddf_b;


    m[0] = ddf_a;
    temp = Chasing(N - 1, Matd, Matc, Mata, d);
    for (int i = 1; i < N; i++)
    {
        m[i] = temp[i - 1];
    }
    m[N] = ddf_b;

    free(d);
    free(Mata);
    free(Matc);
    free(Matd);
    free(h);
    free(temp);
    return m;
}

/**
 * @brief nature Cubic Spline Interpolation.
 * @param N the nubmer of intervals
 * @param a the array on interpolation points
 * @param f the interpolation function
 * @param x the variable where the value of f will return
 * @param ddf_a second order derivative at a
 * @param ddf_b second order derivative at b
 * @return the value of f at x for nature Cubic Spline Interpolation.
 */
double NatureCubicSplineIpl(double (*f)(double), double ddf_a, double ddf_b, double x, int N, double *a)
{
    return CubicSplineIpl(N, f, x, a, NatureCubicSplineIplPara(N, a, f, ddf_a, ddf_b));
}

/**
 * @brief Complete Cubic Spline Interpolation.
 * @param N the nubmer of intervals
 * @param a the array on interpolation points
 * @param f the interpolation function
 * @param x the variable where the value of f will return
 * @param df_a derivative at a
 * @param df_b derivative at b
 * @return the value of f at x for complete Cubic Spline Interpolation.
 */
double CompleteCubicSplineIpl(double (*f)(double), double df_a, double df_b, double x, int N, double *a)
{
    return CubicSplineIpl(N, f, x, a, CompleteCubicSplineIplPara(N, a, f, df_a, df_b));
}

/**
 * @brief Lagrange Cubic Spline Interpolation.
 * @param N the nubmer of intervals
 * @param a the array on interpolation points
 * @param f the interpolation function
 * @param x the variable where the value of f will return
 * @return the value of f at x for Lagrange Cubic Spline Interpolation.
 */
double LagrangeCubicSplineIpl(double (*f)(double), double x, int N, double *a)
{
    return CubicSplineIpl(N, f, x, a, LagrangeCubicSplineIplPara(N, a, f));
}
