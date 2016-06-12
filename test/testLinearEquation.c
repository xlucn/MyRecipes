#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "NR.h"
#include "Test.h"

/**
 * @brief test solving linear equations
 */

static int N = 5;
static double a[] = {0, 1, 1, 1, 1};
static double b[] = {1, -2, 2, -2, 3};
static double c[] = {1, 1, 1, 1, 0};
static double d[] = {2, 4, 4, 4, 4};
static double ans[] = {1, -1, 1, -1, 1};

/**
 * @brief 
 * @returns 
 */
int testChasing()
{
    double *x;

    x = Chasing(N, d, c, a, b);

    if(x == NULL)
    {
        fprintf(stderr, "testChasing: Method failed.");
        return FAILED;
    }
    for(int i = 0; i < N; i++)
    {
        printf("%lf\n", x[i]);
        if(fabs(x[i] - ans[i]) > 1e-8)
        {
            return FAILED;
        }
    }
    return PASSED;
}

static const double LinEqA1[3][3] = {{21, -38, 23},{-38, 70, -43},{23, -43, 27}};
static const double LinEqb1[3] = {-26,49,-28};
static const double LinEqans1[3] = {10.7143,12.5714,9.85714};
static const int LinEqN1 = 3;

/**
 * @brief 
 * @returns 
 */
int testGaussianEli()
{
    double eps = 1e-3;
    
    /* copy */
    int N = LinEqN1;
    double **A = (double **)malloc_s(N * sizeof(double *));
    for(int i = 0; i < N; i++)
    {
        A[i] = (double*)malloc_s(N * sizeof(double));
        for (int j = 0; j < N; j++)
        {
            A[i][j] = LinEqA1[i][j];
        }
    }
    double *b = (double*)malloc_s(N * sizeof(double));
    for(int i = 0; i < N; i++)
    {
        b[i] = LinEqb1[i];
    }


    double *x1 = GaussEli(N, A, b);
    for(int i = 0; i < N; i++)
    {
        printf("%lf\t", x1[i]);
    }
    printf("\n");
    for(int i = 0; i < N; i++)
    {
        if (fabs(x1[i] - LinEqans1[i]) > eps)
        {
            printf("Gauss elimination method failed\n");
            return FAILED;
        }
    }
    return PASSED;
}

/**
 * @brief 
 * @returns 
 */
int testGaussianEliPP()
{
    int N = LinEqN1;
    double eps = 1e-2;
    double **Ae = (double**)malloc_s(N * sizeof(double*));
    for (int i = 0; i < N; i++)
    {
        Ae[i] = (double*)malloc_s((N + 1) * sizeof(double));
        for (int j = 0; j < N; j++)
        {
            Ae[i][j] = LinEqA1[i][j];
        }
        Ae[i][N] = LinEqb1[i];
    }
    double *res = GaussEliPP(N, Ae);
    for(int i = 0; i < N; i++)
    {
        printf("%lf\t", res[i]);
    }
    printf("\n");

    for(int i = 0; i < N; i ++)
    {
        if (fabs(res[i] - LinEqans1[i]) > eps)
        {
            printf("Gauss elimination method with partial pivoting failed\n");
            return FAILED;
        }
    }
    return PASSED;
}

/**
 * @brief 
 * @returns 
 */
int testGaussianEliPPP()
{
    int N = LinEqN1;
    double eps = 1e-2;
    double **Ae = (double**)malloc_s(N * sizeof(double*));
    for (int i = 0; i < N; i++)
    {
        Ae[i] = (double*)malloc_s((N + 1) * sizeof(double));
        for (int j = 0; j < N; j++)
        {
            Ae[i][j] = LinEqA1[i][j];
        }
        Ae[i][N] = LinEqb1[i];
    }
    double *res = GaussEliPPP(N, Ae);
    if(res == NULL)
    {
        return FAILED;
    }
    for(int i = 0; i < N; i++)
    {
        printf("%lf\t", res[i]);
    }
    printf("\n");

    for(int i = 0; i < N; i ++)
    {
        if (fabs(res[i] - LinEqans1[i]) > eps)
        {
            printf("Gauss elimination method with partial pivoting failed\n");
            return FAILED;
        }
    }
    return PASSED;
}
