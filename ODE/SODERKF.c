#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "NR.h"

/**
 * @brief RKF method to solve a system of ODEs
 * @param f a pointer to right-function returning an array of derivatives
 * @param y0 an array of initial values
 * @param a interval
 * @param b 
 * @param m number of equations
 * @param h0 the initial step size
 * @param TOL the required tolerance
 * @param hmax  the minimum and maximum step size allowed
 * @param hmin
 * @param n  temporarily use this number to indicate which method to use
 * @returns struct type SODEsol
 * @note remember to dispose the returned struct by DisposeSODEsol
 * @see DisposeSODEsol
 * @see SODEsol
 */
SODEsol SODERKF(double *(*f)(double, double*), double *y0,
     double a, double b, int m, double h0, double TOL, double hmax, double hmin, int n)
{
    //TODO:Integrate backwards: if a > b then integrate from b to a
    //TODO:make the TOL a array
    if (a > b)
    {
        fprintf(stderr, "the integration region (a,b) not legal\n");
        exit(1);
    }

    long length = 1000; //length of the result, increase by 1000 every time
    long step = 0;
    int TOLflag; // to record if all the residuals are smaller than corresponding TOL
    double h = h0;
    double T = a;
    double *delta = (double *)malloc_s(m * sizeof(double)); // a variable related to the ratio of residuals and TOL
    double *w = (double *)malloc_s(m * sizeof(double)); // a vector for temporary use
    double *R = (double *)malloc_s(m * sizeof(double)); // the residuals
    double **k = (double**)malloc_s(n * sizeof(double*)); // the internal variables
    for(int i = 0; i < n; i++)
    {
        *(k + i) = (double*)malloc_s(m * sizeof(double));
    }
    double *t = (double*)malloc_s(length * sizeof(double)); // the values of variable
    double **y = (double**)malloc_s(length * sizeof(double*)); // the values of the functions
    y[0] = (double*)malloc_s(m * sizeof(double));

    // initialization
    t[0] = T;
    for(int i = 0; i < m; i++)
    {
        w[i] = y0[i];
        y[0][i] = w[i];
    }

    while(T < b)
    {
    	// calculate Ks
        // j, n is for each k variable
        // first we calculate the jth vectork
        for(int j = 0; j < n; j++)
        {
            // the temporary vector parameter w to be passed to functions
            // w = sum( a[j][n]*vectork[n], {n,0,j} )
            for(int icomponent = 0; icomponent < m; icomponent++)
            {
                w[icomponent] = y[step][icomponent];
                for(int indexofks = 0; indexofks < j; indexofks++)
                {
                    w[icomponent] += A78[j * n + indexofks] * k[indexofks][icomponent];
                }
            }
            // i, m is for each component of a variable, the same number as the number of ODEs
            double *temp = f(T + C78[j] * h, w);
            for(int i = 0; i < m; i++)
            {
                // the ith component of vector_k corresponding to the ith function
                // The_jth_vector_k[i] = h*ith_function(t, w)
                k[j][i] = h * temp[i];
            }
        }

        TOLflag = 0;
        // calculate the residual R
        for(int icomponent = 0; icomponent < m; icomponent++)
        {
            R[icomponent] = 0;
            for(int indexofks = 0; indexofks < n; indexofks++)
            {
                R[icomponent] += (B78[indexofks] - Bstar78[indexofks]) * k[indexofks][icomponent] / h;
            }

            delta[icomponent] = pow(TOL / R[icomponent] / 2, 1.0 / 7);

            // lower the step length if the accuracy is not enough
            if(R[icomponent] > TOL)
            {
                h = h * delta[icomponent];
                // printf("The R is %lf, changing step length: %lf\n",R[icomponent], h);
                if(h < hmin)
                {
                    fprintf(stderr, "T = %lf, minimal limit exceeds! lower minimal limit required.\n", T);
                    exit(1);
                }
                TOLflag = 1;
                break;
            }
        }

        // the residuals are all smaller than the TOL, so record the result to array
        if(TOLflag == 0)
        {
            // lengthen the array if the current space is not enough.
            if((step + 1) == length)
            {
                length += 1000;
                t = (double*)realloc_s(t, length * sizeof(double));
                y = (double**)realloc_s(y, length * sizeof(double*));
            }
            // allocate a new line in y.
            y[step+1] = (double*)malloc_s(m * sizeof(double));
            // calculate the new numbers
            for(int icomponent = 0; icomponent < m; icomponent++)
            {
                w[icomponent] = y[step][icomponent];
                for(int indexofks = 0; indexofks < n; indexofks++)
                {
                    w[icomponent] += B78[indexofks] * k[indexofks][icomponent];
                }
                y[step+1][icomponent] = w[icomponent];
            }

            // increase the time t by interval h
            T += h;
            step++;
            t[step] = T;

            double mindelta = 1.5;
            // raise the step length to a suitable value since all the R are below TOL
            for(int icomponent = 0; icomponent < m; icomponent++)
            {
                if(delta[icomponent] < mindelta)
                {
                    mindelta = delta[icomponent];
                }
            }
            if(mindelta > 1.5)
            {
                h = h * mindelta;
            }
        }
    }

    // free the memory allocated before
    free(delta);
    free(w);
    free(R);
    for(int i = 0; i < n; i++)
    {
        free(k[i]);
    }
    free(k);

    // resize the arrays to suitable size -- number of steps
    y = (double**)realloc_s(y, step * sizeof(double*));
    t = (double*)realloc_s(t, step * sizeof(double));
    SODEsol sol = newSODEsol(step, t, y);
    return sol;
}
