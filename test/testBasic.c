#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "NR.h"
#include "Test.h"

/**
 * @brief 
 * @param x 
 * @returns 
 * 
 * 
 */
double f14(double x)
{
    return exp(x);
}

/**
 * @brief 
 * @returns 
 * 
 * 
 */
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


