/**
 * @file testLinearEquation.c
 * @brief test solving linear equations
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "NR.h"
#include "Test.h"

static int N = 5;
static double a[] = {0, 1, 1, 1, 1};
static double b[] = {1, -2, 2, -2, 3};
static double c[] = {1, 1, 1, 1, 0};
static double d[] = {2, 4, 4, 4, 4};
static double ans[] = {1, -1, 1, -1, 1};

/**
 * @brief testChasing
 * @return integer if test passed
 */
int testChasing()
{
    double *x;

    x = Chasing(N, d, c, a, b);

    if(x == NULL)
    {
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

/*--------------------- test general equations -------------------------------*/

/**
 * @brief Linear equation test unit
 */
typedef struct _LinEqtest{
    double *A; /**< empty */
    double *b; /**< empty */
    int N; /**< empty */
    double *ans; /**< empty */
    double eps; /**< empty */
}LinEqtest;

static int _testLinearEquation(double *(*testfunc)(int, double*, double*), LinEqtest t)
{
    double *res = testfunc(t.N, t.A, t.b);
    for(int i = 0; i < t.N; i++)
    {
        printf("%lf\t", res[i]);
        if (fabs(res[i] - t.ans[i]) > t.eps)
        {
            return FAILED;
        }
    }
    printf("\n");
    free(res);
    return PASSED;
}

static double LinEqA1[] = {
    21, -38, 23,
    -38, 70, -43,
    23, -43, 27
};
static double LinEqb1[] = {-26,49,-28};
static double LinEqans1[] = {10.7143,12.5714,9.85714};
static LinEqtest t1 = {LinEqA1, LinEqb1, 3, LinEqans1, 1e-4};

/**
 * @brief testGaussianEli
 * @return integer if test passed
 */
int testGaussianEli()
{
    return _testLinearEquation(GaussEli, t1);
}

/**
 * @brief testGaussianEliPP
 * @return integer if test passed
 */
int testGaussianEliPP()
{
    return _testLinearEquation(GaussEliPP, t1);
}

/**
 * @brief testGaussianEliPPP
 * @return integer if test passed
 */
int testGaussianEliPPP()
{
    return _testLinearEquation(GaussEliPPP, t1);
}

/**
 * @brief testGaussJordanEli
 * @return integer if test passed
 */
int testGaussJordanEli()
{
    return _testLinearEquation(GaussJordanEli, t1);
}

/**
 * @brief testCrout
 * @return integer if test passed
 */
int testCrout()
{
    return _testLinearEquation(Crout, t1);
}
