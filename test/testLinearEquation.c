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



static double LinEqA1[] = {
    21, -38, 23,
    -38, 70, -43,
    23, -43, 27
};
static double LinEqb1[] = {-26,49,-28};
static double LinEqans1[] = {10.7143,12.5714,9.85714};
static int LinEqN1 = 3;

/**
 * @brief 
 * @returns 
 */
int testGaussianEli()
{
    double eps = 1e-3;

    double *x1 = GaussEli(LinEqN1, LinEqA1, LinEqb1);
    for(int i = 0; i < LinEqN1; i++)
    {
        printf("%lf\t", x1[i]);
    }
    printf("\n");
    for(int i = 0; i < LinEqN1; i++)
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
    double eps = 1e-2;
    double *res = GaussEliPP(LinEqN1, LinEqA1, LinEqb1);
    for(int i = 0; i < LinEqN1; i++)
    {
        printf("%lf\t", res[i]);
    }
    printf("\n");

    for(int i = 0; i < LinEqN1; i ++)
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
    double eps = 1e-2;
    double *res = GaussEliPPP(LinEqN1, LinEqA1, LinEqb1);
    if(res == NULL)
    {
        return FAILED;
    }
    for(int i = 0; i < LinEqN1; i++)
    {
        printf("%lf\t", res[i]);
    }
    printf("\n");

    for(int i = 0; i < LinEqN1; i ++)
    {
        if (fabs(res[i] - LinEqans1[i]) > eps)
        {
            printf("Gauss elimination method with partial pivoting failed\n");
            return FAILED;
        }
    }
    return PASSED;
}
