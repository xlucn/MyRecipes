#include <stdlib.h>
#include <math.h>
#include "NR.h"

/**
 * @brief Romberg integration，
 * @param N 所设置的最多序列数,
 * @param [a,b] the interval,
 * @param eps the precision,
 * @param f the integration function
 */
double RombergInt(double(*f)(double), double a, double b, int N, double eps)
{
    double result;
    double h = b - a;
    double **T = (double**)malloc_s(N * sizeof(double*));

    for (int i = 0; i < N; i++)
    {
        *(T + i) = (double*)malloc_s((i + 1) * sizeof(double));
        if (i == 0)
        {
            //T[0][0]
            T[0][0] = (f(b) + f(a)) * h / 2;
        }
        else
        {
            //T[i][0]
            double temp = 0;
            for (int m = 0; m < pow(2, i - 1); m++)
            {
                temp += f(a + (m + 0.5) * h);
            }
            T[i][0] = (T[i - 1][0] + h * temp) / 2;

            //T[i][j]
            for (int j = 1; j < i + 1; j++)
            {
                T[i][j] = (pow(4, j) * T[i][j - 1] - T[i - 1][j - 1]) / (pow(4, j) - 1);
            }

            if (fabs(T[i][i] - T[i][i - 1]) <= eps)
            {
                return T[i][i];
            }

            free(*(T + i - 1));

            h = h / 2;
        }
    }
    result = T[N][N - 1];
    free(*(T + N - 1));
    free(*T);
    return result;

}
