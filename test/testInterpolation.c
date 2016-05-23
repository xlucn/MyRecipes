#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "NR.h"
#include "Test.h"

/**
 * testing interpolation
 */

static double f3(double x)
{
    return 1 / (1 + x * x);
}

/**
 * @brief 
 * @returns 
 */
int testSplineIpl()
{
    double x;
    double a[11] = {-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5};
    double fx[20];
    double SL[20];

    for (int i = 0; i < 20; i++)
    {
        x = 0.5 * (i + 1) - 5.2;
        fx[i] = f3(x);
        SL[i] = NatureCubicSplineIpl(f3, 0.00842, -0.00842,x ,10, a);
        printf("%lf\t%lf\t%lf\n", fx[i], SL[i], SL[i] - fx[i]);
    }
    return PASSED;
}


/**
 * @brief 
 * @returns 
 */
int testHermiteIpl()
{
    int N = 3;
    double a[] = {0, 1, 2};
    double f[] = {1, 2.718, 2.389};
    double df[] = {1, 2.718, 2.389};
    double H;
    double x = 0.25;

    H = Hermite(N, a, f, df, x);

    printf("%lf\n", H);

    return PASSED;
}


