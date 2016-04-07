//LuXu
//Solving linear equations

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "NR.h"
#include "LibFunction.h"

double *GaussEli(int N, double **A, double *b)
{
    double *x = (double *)malloc_s(N * sizeof(double));
    double **l = (double **)malloc_s(N * sizeof(double *));
    for(int i = 0; i < N; i++)
    {
        l[i] = (double *)malloc_s(N * sizeof(double));
    }

    for(int k = 0; k < N - 1; k++)
    {
        for(int i = k + 1; i < N; i++)
        {
            l[i][k] = A[i][k] / A[k][k];
            A[i][k] = 0;
            for(int j = k + 1; j < N; j++)
            {
                A[i][j] -= l[i][k] * A[k][j];
            }
            b[i] -= l[i][k] * b[k];
        }
    }

    for(int k = N - 1; k >= 0; k--)
    {
        x[k] = 0;
        for(int j = k + 1; j < N; j++)
        {
            x[k] += A[k][j] * x[j];
        }
        x[k] = (b[k] - x[k]) / A[k][k];
    }

    for(int i = 0; i < N; i++)
    {
        free(*(l + i));
    }
    free(l);
    return x;
}

/*
Gauss elimination with partial pivoting
*/
double *GaussEliPP(int N, double **a)
{
    int *max = (int*)malloc_s(N * sizeof(int));
    double *temp;
    double *x = (double*)malloc_s(N * sizeof(double));

    for (int k = 0; k < N - 1; k++)
    {
        max[k] = k;
        for (int i = k; i < N; i++)
        {
            if (fabs(a[i][k]) > fabs(a[max[k]][k]))
            {
                max[k] = i;
            }
        }
        if (a[max[k]][k] == 0)
        {
            printf("A is singular\n");
            return NULL;
        }
        if (max[k] != k) {
            temp = a[k];
            a[k] = a[max[k]];
            a[max[k]] = temp;
        }
        for (int i = k + 1; i < N; i++) {
            a[i][k] = a[i][k] / a[k][k];
            for (int j = k + 1; j < N + 1; j++) {
                a[i][j] -= a[i][k] * a[k][j];
            }
        }
    }

    if (a[N - 1][N - 1] == 0) {
        printf("A is singular\n");
        return NULL;
    }
    else
    {
        x[N - 1] = a[N - 1][N] / a[N - 1][N - 1];
    }
    for (int k = N - 2; k >= 0; k--) {
        double t = 0;
        for (int j = k + 1; j < N; j++) {
            t += a[k][j] * x[j];
        }
        x[k] = (a[k][N] - t) / a[k][k];
    }
    free(max);
    return x;
}

/*
Gauss elimination with partial pivoting proportionally
*/
double *GaussEliPPP(int N, double **a)
{
    int r = 0;
    double q, *t;
    double *max = (double*)malloc_s(N * sizeof(double));
    double *x = (double*)malloc_s(N * sizeof(double));

    // get the largest number of each line
    for (int i = 0; i < N; i++)
    {
        max[i] = 0;
        for (int j = 0; j < N; j++)
        {
            max[i] = (fabs(max[i]) < fabs(a[i][j])) ? fabs(a[i][j]) : max[i];
        }
        if (max[i] == 0)
        {
            printf("A is singular\n");
        }
    }

    // do the gauss elimination with partial pivoting proportionally
    for (int k = 0; k < N - 1; k++)
    {
        // choose the pivot element
        for (int i = k; i < N; i++)
        {
            r = (fabs(a[r][k] / max[r]) < fabs(a[i][k] / max[i])) ? i : r;
        }
        if (a[r][k] == 0)
        {
            printf("A is singular\n");
            return NULL;
        }
        if (r != k)
        {
            q = max[k]; max[k] = max[r]; max[r] = q;
            t = a[k]; a[k] = a[r]; a[r] = t;
        }

        for (int i = k + 1; i < N; i++)
        {
            a[i][k] = a[i][k] / a[k][k];
            for (int j = k + 1; j < N + 1; j++)
            {
                a[i][j] -= a[i][k] * a[k][j];
            }
        }
    }

    if (a[N - 1][N - 1] == 0)
    {
        printf("A is singular\n");
        return NULL;
    }

    for (int k = N - 1; k >= 0; k--)
    {
        x[k] = 0;
        for (int j = k + 1; j < N; j++)
        {
            x[k] += a[k][j] * x[j];
        }
        x[k] = (a[k][N] - x[k]) / a[k][k];
    }
    free(max);
    return x;
}

double *GaussJordanEli(int N, double **a)
{
    int* max = (int*)malloc_s(N * sizeof(int));
    double* x = (double*)malloc_s(N * sizeof(double));

    for (int k = 0; k < N; k++)
    {
        max[k] = 0;
        for (int i = 0; i < N; i++)
        {
            int flag = 1;
            for (int j = 0; j < k; j++)
            {
                if (i == max[j])
                {
                    flag = 0;
                }
            }
            if (flag)
            {
                max[k] = (fabs(a[max[k]][k]) < fabs(a[i][k])) ? i : max[k];
            }
        }
        if (a[max[k]][k] == 0)
        {
            printf("A is singular\n");
            return NULL;
        }

        for (int i = 0; i < N; i++)
        {
            if (i != max[k])
            {
                a[i][k] = a[i][k] / a[max[k]][k];
                for (int j = k + 1; j < N + 1; j++)
                {
                    a[i][j] -= a[i][k] * a[max[k]][j];
                }
            }

        }
    }

    for (int k = 0; k < N; k++)
    {
        x[k] = a[max[k]][N] / a[max[k]][k];
    }
    free(max);
    return x;
}

double *Chasing(int N, double *d, double *c, double *a, double *b)
{
    double *p = (double *)malloc_s(N * sizeof(double));
    double *q = (double *)malloc_s(N * sizeof(double));
    double *x = (double *)malloc_s(N * sizeof(double));
    double *y = (double *)malloc_s(N * sizeof(double));

    if(d[0] == 0)
    {
        return NULL;
    }
    p[0] = d[0];
    q[0] = c[0] / d[0];
    for(int i = 1; i < N - 1; i++)
    {
        p[i] = d[i] - a[i] * q[i - 1];
        if(p[i] == 0)
        {
            return NULL;
        }
        q[i] = c[i] / p[i];
    }
    p[N - 1] = d[N - 1] - a[N - 1] * q[N - 2];
    if(p[N - 1] == 0)
    {
        return NULL;
    }
    y[0] = b[0] / p[0];
    for(int i = 1; i < N; i++)
    {
        y[i] = (b[i] - a[i] * y[i - 1]) / p[i];
    }
    x[N - 1] = y[N - 1];
    for(int i = N - 2; i >= 0; i--)
    {
        x[i] = y[i] - q[i] * x[i + 1];
    }
    free(p);
    free(q);
    free(y);
    return x;
}
