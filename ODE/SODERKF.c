/** @file SODERKF.c */
#include <stdio.h>
#include <math.h>
#include "NR.h"
#include "constants.h"
#include "NRprivate.h"

static SODEsol SODERKF(double *(*f)(double, double*), double *y0,
     double a, double b, int m, double h, double TOL, double hmax, double hmin, int n)
{
    const double *A, *B, *Bstar, *C;
    if(n == 13)
        A = A78; B = B78; Bstar = Bstar78; C = C78;
    else if(n == 6)
        A = A45; B = B45; Bstar = Bstar45; C = C45;

    if (a > b)
    {
        fprintf(stderr, "the integration region (a,b) not legal\n");
        exit(1);
    }
    int length = 1024;
    int step = 0;
    double *delta = newArray1d(m);
    double *w = newArray1d(m); // a vector for temporary use
    double *R = newArray1d(m); // the residuals
    double **k = (double**)malloc_s(n * sizeof(double*));

    double *t = newArray1d(length);
    double **y = (double**)malloc_s(length * sizeof(double*));
    /* initialization */
    t[0] = a;
    y[0] = newArray1d(m);
    for(int i = 0; i < m; i++)
    {
        y[0][i] = y0[i];
    }

    while(t[step] < b)
    {
    	/* calculate Ks */
        for(int in = 0; in < n; in++)
        {
            for(int im = 0; im < m; im++)
            {
                w[im] = y[step][im];
                for(int j = 0; j < in; j++)
                {
                    w[im] += A[in * n + j] * k[j][im] * h;
                }
            }
            k[in] = f(t[step] + C[in] * h, w);
        }

        int TOLflag = 0; /* to record if all the residuals are smaller than corresponding TOL */
        /* calculate the residual R */
        for(int im = 0; im < m; im++)
        {
            R[im] = 0;
            for(int in = 0; in < n; in++)
            {
                R[im] += (B[in] - Bstar[in]) * k[in][im];
            }
            delta[im] = pow(TOL / fabs(R[im]) / 2, 1.0 / 7);

            /* lower the step length if the accuracy is not enough */
            if(R[im] > TOL)
            {
                h *= delta[im] <= 0.1 ? 0.1 : (delta[im] >= 1.0 ? 1.0 : delta[im]);
                if(h < hmin)
                {
                    fprintf(stderr, "T = %lf, minimal limit exceeds! lower minimal limit required.\n", t[step]);
                    exit(1);
                }
                TOLflag = 1;
                break;
            }
        }

        /* the residuals are all smaller than the TOL, so record the result to array */
        if(TOLflag == 0)
        {
            step++;
            /* lengthen the array if the current space is not enough. */
            if(step == length)
            {
                length += 1024;
                t = (double*)realloc_s(t, length * sizeof(double));
                y = (double**)realloc_s(y, length * sizeof(double*));
            }
            y[step] = newArray1d(m);
            /* calculate the new numbers */
            for(int im = 0; im < m; im++)
            {
                y[step][im] = y[step - 1][im];
                for(int in = 0; in < n; in++)
                {
                    y[step][im] += B[in] * k[in][im] * h;
                }
            }
            t[step] = t[step - 1] + h;

            double mindelta = 1;
            for(int im = 0; im < m; im++) if(delta[im] < mindelta)
            {
                mindelta = delta[im];
            }
            h *= mindelta <= 1.0 ? 1.0 : (mindelta >= 2.0 ? 2.0 : mindelta);
            h = h > hmax ? hmax : h;
        }
    }

    // free the memory allocated before
    delArray1d(delta);
    delArray1d(w);
    delArray1d(R);
    free(k);

    // resize the arrays to suitable size -- number of steps
    y = (double**)realloc_s(y, step * sizeof(double*));
    t = Array1dResize(t, step);
    SODEsol sol = newSODEsol(step, t, y);
    return sol;
}
/**
 * @brief RKF method to solve a system of ODEs
 * @param f a pointer to right-function returning an array of derivatives
 * @param y0 an array of initial values
 * @param a lower limit of interval
 * @param b upper limit of interval
 * @param m number of equations
 * @param h the initial step size
 * @param TOL the required tolerance
 * @param hmax the minimum step size allowed
 * @param hmin the maximum step size allowed
 * @returns struct type SODEsol
 * @note remember to dispose the returned struct by DisposeSODEsol
 * @see DisposeSODEsol
 * @see SODEsol
 * @todo Integrate backwards: if a > b then integrate from b to a
 * @todo make the TOL a array
 */
SODEsol SODERKF45(double *(*f)(double, double*), double *y0,
     double a, double b, int m, double h, double TOL, double hmax, double hmin)
{
    return SODERKF(f, y0, a, b, m, h, TOL, hmax, hmin, 6);
}
/**
 * @brief RKF method to solve a system of ODEs
 * @param f a pointer to right-function returning an array of derivatives
 * @param y0 an array of initial values
 * @param a lower limit of interval
 * @param b upper limit of interval
 * @param m number of equations
 * @param h the initial step size
 * @param TOL the required tolerance
 * @param hmax the minimum step size allowed
 * @param hmin the maximum step size allowed
 * @returns struct type SODEsol
 * @note remember to dispose the returned struct by DisposeSODEsol
 * @see DisposeSODEsol
 * @see SODEsol
 * @todo Integrate backwards: if a > b then integrate from b to a
 * @todo make the TOL a array
 */
SODEsol SODERKF78(double *(*f)(double, double*), double *y0,
     double a, double b, int m, double h, double TOL, double hmax, double hmin)
{
    return SODERKF(f, y0, a, b, m, h, TOL, hmax, hmin, 13);
}
