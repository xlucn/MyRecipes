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

static double LinEqA1[3][3] = {{21, -38, 23},{-38, 70, -43},{23, -43, 27}};
static double LinEqb1[3] = {-26,49,-28};
static double LinEqans[] = {10.7143,12.5714,9.85714};
static int LinEqN = 3;

/**
 * @brief 
 * @returns 
 */
int testGaussianEli()
{
    int N = LinEqN;
    double eps = 1e-3;
    double **A = (double **)malloc_s(N * sizeof(double *));
    for(int i = 0; i < N; i++)
    {
        A[i] = (double*)malloc_s(N * sizeof(double));
        for (int j = 0; j < N; j++)
        {
            A[i][j] = LinEqA1[i][j];
        }
    }


    double *x1 = GaussEli(N, A, LinEqb1);
    for(int i = 0; i < N; i++)
    {
        printf("%lf\t", x1[i]);
    }
    printf("\n");
    for(int i = 0; i < N; i++)
    {
        if (fabs(x1[i] - LinEqans[i]) > eps)
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
    int N = LinEqN;
    double eps = 1e-4;
    double **Ae = (double**)malloc_s(N * sizeof(double*));
    for (int i = 0; i < N; i++)
    {
        Ae[i] = (double*)malloc_s((N + 1) * sizeof(double));
        for (int j = 0; j < N; j++)
        {
            Ae[i][j] = LinEqA1[i][j];
        }
        Ae[i][N] = b[i];
    }
    double *x2 = GaussEliPP(LinEqN, Ae);
    for(int i = 0; i < N; i++)
    {
        printf("%lf\t", x2[i]);
    }
    printf("\n");

    for(int i = 0; i < N; i ++)
    {
        if (fabs(x2[i] - LinEqans[i]) > eps)
        {
            printf("Gauss elimination method with partial pivoting failed\n");
            return FAILED;
        }
    }
    return PASSED;
}

