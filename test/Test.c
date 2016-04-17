/**
 * LuXu
 * test functions
 * functions needed by each testfunction is defined right before it.
 *
 * note: if a function is declared as "int test...()", it will be collected by
 *       GenerateTest.py and regarded as a test function.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "NR.h"
#include "LibFunction.h"
#include "Test.h"

void testall()
{
    // the num is defined in "Test.h", which is the number of test functions.
    int* result = (int*)malloc_s(num * sizeof(int));  // record the test results

    // the tests is defined in "Test.h", which is the list of test functions.
    for(int i = 0; i < num; i ++)
    {
        printf("=========================================\n");
        printf("* * * * * * * *testing No.%d* * * * * * *\n", i);
        result[i] = tests[i]();
        printf("=========================================\n\n");
    }

    for (int i = 0; i < num; i++)
    {
        printf("%s No.%3d %s\n", (result[i]) ? "*[failed]" : " [passed]", i, names[i]);
    }
}

/**
 * testing integration
 */
static double integrand1(double x)
{
    return x * sqrt(1 + x * x);
}

static double integral1(double x)
{
    return  pow(1 + x * x, 1.5) / 3;
}

static double integrand2(double x)
{
    return 1 / (1 + x);
}

static double integral2(double x)
{
    return log(x + 1);
}

int testAdaptiveSimpson()
{
    double a = 0;
    double b = 3;
    double TOL = 1e-4;
    double res = AdaptiveSimpsonInt(integrand1, a, b, TOL);
    double ans = integral1(b) - integral1(a);
    printf("T\t\t = %.12lf\nI\t\t = %.12lf\ndelta\t = %.12lf\n",
        res, ans, fabs(res - ans));

    return fabs(res - ans) < TOL ? PASSED : FAILED;
}

int testRomberg()
{
    int N = 100;
    double a = 0;
    double b = 1.5;
    double TOL = 1e-10;
    double res = RombergInt(integrand2, a, b, N, TOL);
    double ans = integral2(b) - integral2(a);
    printf("T\t = %.12lf\nI\t = %.12lf\ndelta\t = %.12lf\n",
        res, ans, fabs(res - ans));

    return fabs(res - ans) < TOL ? PASSED : FAILED;
}

/**
 * testing solving linear equations
 */

static int N = 5;
static double a[] = {0, 1, 1, 1, 1};
static double b[] = {1, -2, 2, -2, 3};
static double c[] = {1, 1, 1, 1, 0};
static double d[] = {2, 4, 4, 4, 4};
static double ans[] = {1, -1, 1, -1, 1};

int testChasing()
{
    double *x;

    x = Chasing(N, d, c, a, b);

    if(x == NULL)
    {
        fprintf(stderr, "testChasing: Method failed.");
        return FAILED;
    }
    else
    {
        for(int i = 0; i < N; i++)
        {
            printf("%f\n", x[i]);
            if(fabs(x[i] - ans[i]) > 1e-8)
            {
                return FAILED;
            }
        }
    }
    return PASSED;
}

static double LinEqA1[3][3] = {{21, -38, 23},{-38, 70, -43},{23, -43, 27}};
static double LinEqb1[3] = {-26,49,-28};
static double LinEqans[] = {10.7143,12.5714,9.85714};
static int LinEqN = 3;

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


/**
 * testing solving ODEs
 */

static double f4(double t, double y)
{
    return (y * y + y) / t;
}

static double y4(double t)
{
    return 2 * t / (1 - 2 * t);
}

int testClassicRK()
{
    double t;
    int N = 4;
    double y0 = -2;
    double a = 1;
    double b = 3;

    double *result = ClassicRungeKutta(f4, a, b, y0, N);
    for(int i = 0; i < N + 1; i ++)
    {
        t = a + i * (b - a) / N;
        printf("%lf\t%lf\n", result[i], y4(t));
    }
    return PASSED;
}


static double f5(double t, double y)
{
    return -y + t + 1;
}

static double y5(double t)
{
    return exp(-t) + t;
}

int testAdamsPECE()
{
    double a = 0;
    double b = 1;
    int N = 5;
    double y0 = y5(a);
    double dy0 = f5(a, y0);
    double *ans = (double*)malloc_s((N + 1) * sizeof(double));
    for(int i = 0; i < N + 1; i++)
    {
        ans[i] = y5(a + (b - a) * i / N);
    }

    double *result = AdamsPECE(f5, a, b, dy0, y0, N);
    for(int i = 0; i < N + 1; i ++)
    {
        printf("%lf\t%lf\n", result[i], ans[i]);
        if(fabs(result[i] - ans[i]) > 0.01)
        {
            return FAILED;
        }
    }
    return PASSED;
}

static double f6(double t, double y)
{
    return - y + t * t + 1;
}

static double y6(double t)
{
    return - 2 / exp(t) + 3 - 2 * t + t * t;
}

int testRKF()
{
    double a = 0;
    double b = 1;
    double TOL = 1e-5;
    double hmax = 1e-2;
    double hmin = 1e-6;
    double y0 = 1;
    double *result;

    result = RKF78(f6, a, b, y0, TOL, hmax, hmin);

    printf("Testing RKF Method\n");
    printf("t\t\th\t\ty\n");
    for(int i = 0; i < result[0]; i ++)
    {
        for(int j = 0; j < 3; j++)
        {
            printf("%lf\t", result[i * 3 + j + 1]);
        }
        printf("\n");

        // Check the results
        if (fabs(result[i * 3 + 3] - y6(result[i * 3 + 1])) > 1e-5)
        {
            return FAILED;
        }
    }
    return PASSED;
}

static double SODEf1(double t, double *y)
{
    return -4 * y[0] - 2 * y[1] + cos(t) + 4 * sin(t);
}

static double SODEf2(double t, double *y)
{
    return 3 * y[0] + y[1] - 3 * sin(t);
}

static double SODEy1(double x)
{
    return exp(-2 * x) * (-2 + 2 * exp(x) + exp(2 * x) * sin(x));
}

static double SODEy2(double x)
{
    return -exp(-2 * x) * (3 * exp(x) - 2);
}


int testSODERungeKutta()
{
    double (*f[2])(double,double*) = {SODEf1, SODEf2};
    double (*y[2])(double) = {SODEy1, SODEy2};
    double a = 0;
    double b = 1;
    double y0[] = {0, -1};
    int N = 10;
    int m = 2;
    double **result = SODERungeKutta(f, a, b, y0, m, N);

    double *t;
    double **res;
    SODERKF(&t, &res, f, y0, a, b, m, 0.1, 1e-8, 1e4, 1e-4, 13);

    printf("Testing Runge-Kutta Method for a System of ODEs\n");
    printf("%8s%16s%16s%16s%16s\n", "t", "result1", "y1", "result2", "y2");
    for(int i = 0; i < N + 1; i ++)
    {
        printf("%8.2lf", a + i * (b - a) / N);
        for(int j = 0 ; j < m; j++)
        {
            printf("%16lf%16lf", result[i][j], y[j](a + i * (b - a) / N));
        }
        printf("\n");
    }
    return PASSED;
}


/**
 * testing orthogonal decomposition
 */

int testLeastSq()
{
    int m = 4;
    int n = 3;
    double temp[4][3] = {{1, -2, 1},{0, 1, -1},{2, -4, 3},{4, -7, 4}};
    double **A = (double**)malloc_s(m * sizeof(double));
    for (int i = 0; i < m; i++)
    {
        A[i] = temp[i];
    }
    double b[4] = {-4, 3, 1, -6};

    double *res = LeastSquare(m, n, A, b);
    for (int i = 0; i < n; i++)
    {
        printf("%lf\t", res[i]);
    }
    printf("\n");
    return PASSED;
}

int testQR()
{
    int m = 4;
    int n = 3;
    double eps = 1e-6;
    double temp[4][3] = {{1,-1,1},{0,-2,-1},{1,1,1},{0,2,1}};
    double **A = (double**)malloc_s(m * sizeof(double*));
    for(int i = 0; i < m; i++)
    {
        A[i] = temp[i];
    }

    double ***QR = ImprovedGramSchmidtQR(m, n, A);
    for(int i = 0; i < m; i++)
    {
        for(int j = 0; j < n; j++)
        {
            printf("%lf\t",QR[0][i][j]);
        }
        printf("\n");
    }
    printf("\n");
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            printf("%lf\t",QR[1][i][j]);
        }
        printf("\n");
    }
    printf("\n");
    for(int i = 0; i < m; i ++)
    {
        for(int j = 0; j < n; j++)
        {
            for(int k = 0; k < n; k++)
            {
                A[i][j] -= QR[0][i][k] * QR[1][k][j];
            }
            if(fabs(A[i][j]) > eps)
            {
                return FAILED;
            }
        }
    }
    return PASSED;
}


/**
 * testing solving a equations
 */

static double f9(double x)
{
    return x * x * x - x - 1;
}

int testBisection()
{
    double a = 1;
    double b = 2;
    double eps = 1e-8;
    double p;

    p = Bisection(f9, a, b, eps);
    printf("%lf\n", p);
    if(fabs(p - 1.32471795724475) > eps)
    {
        return FAILED;
    }
    return PASSED;
}

static double f10(double x)
{
    return (x + 2) * x * x - 4;
}
static double g10(double x)
{
    return x - f10(x) / (3 * x + 4) / x;
}

int testPicardRecurtion()
{
    double x = 1;
    double eps = 1e-4;
    double p;

    p = PicardIteration(g10, x, eps);
    printf("%lf\n", p);
    if(fabs(p - 1.130395435) > eps)
    {
        return FAILED;
    }
    return PASSED;
}

static double g11(double x)
{
    return sqrt(4 / (2 + x));
}

int testSteffensen()
{
    double x = 1.5;
    double eps = 1e-8;
    double p;

    p = SteffensenIteration(g11, x, eps);
    printf("%lf\n", p);
    if(fabs(p - 1.130395435) > eps)
    {
        return FAILED;
    }
    return PASSED;
}

static double f12(double x)
{
    return ((x + 2) * x + 10) * x - 20;
}

static double df12(double x)
{
    return (3 * x + 4) * x + 10;
}

int testNewtonMethod()
{
    double x = 1;
    double eps = 1e-8;
    double p;

    p = NewtonMethod(f12, df12, x, eps);
    printf("%lf\n", p);
    if(fabs(p - 1.368808108) > eps)
    {
        return FAILED;
    }
    return PASSED;
}

double f13(double x)
{
    return (2 * x * x - 5) * x - 1;
}

int testSecent()
{
    double x0 = 2;
    double x1 = 1;
    double eps = 1e-8;
    double p;

    p = SecentMethod(f13, x0, x1, eps);
    printf("%lf\n", p);
    if(fabs(p - 1.67298165) > eps)
    {
        return FAILED;
    }
    return PASSED;
}

double f14(double x)
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

int testMuller()
{
    return PASSED;
}
