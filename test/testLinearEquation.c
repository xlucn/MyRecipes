#include <stdio.h>
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

typedef struct _LinEqtest{
    double *A;
    double *b;
    int N;
    double *ans;
    double eps;
}LinEqtest;

/**
 * @brief 
 * @returns 
 */
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
    delArray1d(res);
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
 * @brief 
 * @returns 
 */
int testGaussianEli()
{
    return _testLinearEquation(GaussEli, t1);
}

/**
 * @brief 
 * @returns 
 */
int testGaussianEliPP()
{
    return _testLinearEquation(GaussEliPP, t1);
}

/**
 * @brief 
 * @returns 
 */
int testGaussianEliPPP()
{
    return _testLinearEquation(GaussEliPPP, t1);
}

int testCrout()
{
    return _testLinearEquation(Crout, t1);
}
