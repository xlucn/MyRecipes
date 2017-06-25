/** @file testInterpolation.c */
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



static double f14(double x)
{
    return exp(x);
}

int testDividedDiff()
{
    int N = 7;
    int k = 4;
    double eps = 1e-3;
    double x[] = {1, 2, 3, 4, 5, 6, 7};
    double **d = FullDividedDiff(f14, x, N, k);
    double **m = DividedDiffMatrix(f14, x, N);
    double answer[7][7] =
    {
        {2.718282,	7.389056,	20.085537,	54.598150,	148.413159,	403.428793,	1096.633158},
        {4.670774,	12.696481,	34.512613,	93.815009,	255.015634,	693.204365},
        {4.012853,	10.908066,	29.651198,	80.600313,	219.094365},
        {2.298404,	6.247711,	16.983038,	46.164684},
        {0.987327,	2.683832,	7.295411},
        {0.339301,	0.922316},
        {0.097169}
    };


    for (int i = 0; i < k; i++)
    {
        for (int j = 0; j < N - i; j++)
        {
            if (fabs(d[i][j] - answer[i][j]) > eps) {
                printf("failed\n");
                return FAILED;
            }
        }
    }

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N - i; j++)
        {
            if (fabs(m[j][i + j] - answer[i][j]) > eps) {
                printf("failed\n");
                return FAILED;
            }
        }
    }
    return PASSED;
}
