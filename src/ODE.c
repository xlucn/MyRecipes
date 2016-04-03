//2015.11.27
//卢旭
//Solving Ordinary differential equations
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "NumericalRecipes.h"
#include "LibFunction.h"

/*
Euler Method to solve ODE. y0:initial value, f(t, y) = dy/dt
*/
double* Euler(double f(double, double), double a, double b, double y0, int N)
{
    double h = (b - a) / N;
    double x = a;
    double* y = (double*)malloc_s((N + 1) * sizeof(double));
    y[0] = y0;

    for(int i = 0; i < N; i++)
    {
        y[i + 1] += h * f(x, y[i]);
        x += h;
    }
    return y;
}

/*
2阶方法的一般步骤
*/
static double* TwoStageRungeKutta(int N, double y0, double a, double b, double f(double, double), double para)
{
    double c1 = 1 -  0.5 / para;
    double c2 = 1 - c1;
    double k1;
    double k2;
    double h = (b - a) / N;
    double x = a;
    double y = y0;
    double* result = (double*)malloc_s(N * sizeof(double));

    for(int i = 0; i < N; i++)
    {
        k1 = f(x, y);
        k2 = f(x + para * h, y + para * h * k1);
        y += h * (c1 * k1 + c2 * k2);
        x += h;
        result[i] = y;
    }
    return result;
}
/*
改善的欧拉方法
*/
double* ImprovedEuler(double f(double, double), double a, double b, double y0, int N)
{
    return TwoStageRungeKutta(N, y0, a, b, f, 1);
}
/*
中点方法或变形的欧拉方法
*/
double* MID(double f(double, double), double a, double b, double y0, int N)
{
    return TwoStageRungeKutta(N, y0, a, b, f, 1 / 2);
}
/*
Heun方法
*/
double* Heun(double f(double, double), double a, double b, double y0, int N)
{
    return TwoStageRungeKutta(N, y0, a, b, f, 2 / 3);
}

/*
三阶Heun方法
*/
double *ThreeStageHeun(int N, double y0, double a, double b, double f(double, double))
{
    return 0;
}

/*
ThreeStageRungeKuttaMathod
*/
double *ThreeStageRungeKuttaMathod(double f(double, double), double a, double b, double y0, int N)
{
    return 0;
}

/*
Classic Runge-Kutta Method
*/
double *ClassicRungeKutta(double f(double, double), double a, double b, double y0, int N)
{
    double k[4];
    double h = (b - a) / N;
    double x = a;
    double y = y0;
    double* result = (double*)malloc_s((N + 1) * sizeof(double));
    result[0] = y0;

    for(int i = 0; i < N; i++)
    {
        k[0] = f(x, y);
        k[1] = f(x + h / 2, y + h * k[0] / 2);
        k[2] = f(x + h / 2, y + h * k[1] / 2);
        k[3] = f(x + h, y + h * k[2]);
        y += h / 6 * (k[0] + 2 * k[1] + 2 * k[2] + k[3]);
        x += h;
        result[i + 1] = y;
    }

    return result;
}

static double A78[13][13] = {
    {0},
    {2.0/17.0, 0},
    {1.0/36.0, 1.0/12.0, 0},
    {1.0/24.0, 0, 1.0/8.0, 0},
    {5.0/12.0, 0, -25.0/16.0, 25.0/16.0, 0},
    {1.0/20.0, 0, 0, 1.0/4.0, 1.0/5.0, 0},
    {-25.0/108.0, 0, 0, 125.0/108.0, -65.0/27.0, 125.0/54.0, 0},
    {31.0/300.0, 0, 0, 0, 61.0/225.0, -2.0/9.0, 13.0/900.0, 0},
    {2.0, 0, 0, -53.0/6.0, 704.0/45.0, -107.0/9.0, 67.0/90.0, 3.0, 0},
    {-91.0/108.0, 0, 23.0/108.0, -976.0/135.0, 311.0/54.0, -19.0/60.0, 17.0/6.0, -1.0/12.0, 0},
    {2383.0/4100.0, 0, 0, -341.0/164.0, 4496.0/1025.0, -301.0/82.0, 2133.0/4100.0, 45.0/82.0, 45.0/164.0, 18.0/41.0, 0},
    {3.0/205.0, 0, 0, 0, 0, -6.0/41.0, -3.0/205.0, -3.0/41.0, 3.0/41.0, 6.0/41.0, 0, 0},
    {-1777.0/4100.0, 0, 0, -341.0/164.0, 4496.0/1025.0, -289.0/82.0, 2193.0/4100.0, 51.0/82.0, 33.0/164.0, 12.0/41.0, 0, 1.0, 0}
};
static double B78[13] = {41.0/840.0, 0, 0, 0, 0, 34.0/105.0, 9.0/35.0, 9.0/35.0, 9.0/280.0, 9.0/280.0, 41.0/840.0, 0, 0};
static double Bstar78[13] = {0, 0, 0, 0, 0, 34.0/105.0, 9.0/35.0, 9.0/35.0, 9.0/280.0, 9.0/280.0, 0, 41.0/840.0, 41.0/840.0};
static double C78[13] = {0, 2.0/27.0, 1.0/9.0, 1.0/6.0, 5.0/12.0, 1.0/2.0, 5.0/6.0, 1.0/6.0, 2.0/3.0, 1.0/3.0, 1.0, 0, 1.0};

static double A45[6][6] =
{
    {0},
    {1.0/4.0, 0},
    {3.0/32.0, 9.0/32.0, 0},
    {1932.0/2197.0, -7200.0/2197.0, 7296.0/2197.0, 0},
    {439.0/216.0, -8.0, 3680.0/513.0, -845.0/4104.0, 0},
    {-8.0/27.0, 2.0, -3544.0/2565.0, 1859.0/4104.0, -11.0/40.0, 0}
};
static double B45[6] = {25.0/216.0, 0, 1408.0/2565.0, 2197.0/4104.0, -1.0/5.0, 0};
static double Bstar45[6] = {16.0/135.0, 0, 6656.0/12825.0, 28561.0/56430.0, -9.0/50.0, 2.0/55.0};
static double C45[6] = {0, 1.0/4.0, 3.0/8.0, 12.0/13.0, 1, 1.0/2.0};


static double *RKFmn(double f(double,double), double a, double b, double y0, double TOL, double hmax, double hmin,
    double** A, double* B, double* Bstar, double* C, int n)
{
    int step = 0; //the total steps
    double t = a;
    double y = y0;
    double h = hmax;
    double *k = (double*)malloc_s(n * sizeof(double));
    double R = 0;
    double delta;
    double *result;
    step++;
    result = (double*)malloc_s((step * 3 + 1) * sizeof(double));
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
                y += A[i][j] * k[j];
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

double *RKF78(double f(double,double), double a, double b, double y0, double TOL, double hmax, double hmin)
{
    double **A = (double**)malloc_s(13 * sizeof(double*));
    for(int i = 0; i < 13; i++)
    {
        A[i]=A78[i];
    }
    return RKFmn(f, a, b, y0, TOL, hmax, hmin, A, B78, Bstar78, C78, 13);
}

double *RKF45(double f(double,double), double a, double b, double y0, double TOL, double hmax, double hmin)
{
    double **A = (double**)malloc_s(6 * sizeof(double*));
    for(int i = 0; i < 6; i++)
    {
        A[i]=A45[i];
    }
    return RKFmn(f, a, b, y0, TOL, hmax, hmin, A, B45, Bstar45, C45, 6);
}


/*
Adams显式和隐式方法的PECE模式校正方法，这里k=1，用经典Runge-Kutta方法提供初值
*/
double *AdamsPECE(double f(double, double), double a, double b, double dy0, double y0, int N)
{
    double *y = (double *)malloc_s((N + 1) * sizeof(double));
    double *dy = (double *)malloc_s((N + 1) * sizeof(double));
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
    free(dy);
    return y;
}

double **SODERungeKutta(double (**f)(double, double*), double a, double b, double *y0, int m, int N)
{
    double h = (b - a) / N;
    double t = a;
    double *w = (double *)malloc_s(m * sizeof(double));
    double **y = (double**)malloc_s((N + 1) * sizeof(double*));
    for(int i = 0; i < N + 1; i++)
    {
        *(y + i) = (double*)malloc_s(m * sizeof(double));
    }
    double **k = (double**)malloc_s(4 * sizeof(double*));
    for(int i = 0; i < 4; i++)
    {
        *(k + i) = (double*)malloc_s(m * sizeof(double));
    }

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
    return y;
}


// m: number of functions or ODEs,
// n: number of ks, order of RKF method matrix
int SODERKF(double **t, double ***y, double (**f)(double, double*), double *y0,
    double a, double b, int m, double h0, double TOL, double hmax, double hmin, int n)
{
    long length = 1000;
    long step = 0;
    int TOLflag; // to record if all the residuals are smaller than recosponding TOL
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
    *t = (double*)malloc_s(length * sizeof(double)); // the values of variable
    *y = (double**)malloc_s(length * sizeof(double*)); // the values of the functions
    (*y)[0] = (double*)malloc_s(m * sizeof(double));

    // initialization
    (*t)[0] = T;
    for(int i = 0; i < m; i++)
    {
        w[i] = y0[i];
        (*y)[0][i] = w[i];
    }

    while(T < b)
    {
        // j, n is for each k variable
        // first we calculate the jth vectork
        for(int j = 0; j < n; j++)
        {
            // the temporary vector parameter w to be passed to functions
            // w = sum( a[j][n]*vectork[n], {n,0,j} )
            for(int icomponent = 0; icomponent < m; icomponent++)
            {
                w[icomponent] = (*y)[step][icomponent];
                for(int indexofks = 0; indexofks < j; indexofks++)
                {
                    w[icomponent] += A78[j][indexofks] * k[indexofks][icomponent];
                }
            }
            // i, m is for each component of a variable, the same number as the number of ODEs
            for(int i = 0; i < m; i++)
            {
                // the ith component of vectork corresponding to the ith function
                // The_jth_vectork[i] = h*ith_function(t, w)
                k[j][i] = h * f[i](T + C78[j] * h, w);
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
                // printf("The R is %f, changing step length: %f\n",R[icomponent], h);
                if(h < hmin)
                {
                    printf("minimal limit exceeds! lower minimal limit required.\n");
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
                (*t) = (double*)realloc_s(*t, length * sizeof(double));
                (*y) = (double**)realloc_s(*y, length * sizeof(double*));
            }
            // allocate a new line in y.
            (*y)[step+1] = (double*)malloc_s(m * sizeof(double));
            // use the temporary vector w to temporarily assign the result
            for(int icomponent = 0; icomponent < m; icomponent++)
            {
                w[icomponent] = (*y)[step][icomponent];
                for(int indexofks = 0; indexofks < n; indexofks++)
                {
                    w[icomponent] += B78[indexofks] * k[indexofks][icomponent];
                }
                (*y)[step+1][icomponent] = w[icomponent];
            }

            // increase the time t by interval h
            T += h;
            step++;
            (*t)[step] = T;

            double mindelta = 1;
            // raise the step length to a suitable value since all the R are below TOL
            for(int icomponent = 0; icomponent < m; icomponent++)
            {
                if(delta[icomponent] < mindelta)
                {
                    mindelta = delta[icomponent];
                }
            }
            if(mindelta > 1)
            {
                h = h * mindelta;
                // printf("The delta is %f, changing step length: %f\n",mindelta, h);
            }
        }
    }


    *y = (double**)realloc_s(*y, step * sizeof(double*));
    *t = (double*)realloc_s(*t, step * sizeof(double));
    return step;
}
