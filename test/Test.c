/**
 * @file Test.c
 * @author Lu Xu
 * @brief test functions
 * @note functions needed by each test function is defined right before it.
 *
 * @note if a function is declared as "int test...()", it will be collected by
 *       GenerateTest.py and regarded as a test function.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "NR.h"
#include "LibFunction.h"
#include "Test.h"

/**
 * @brief Test all the functions in this project
 * 
 * This testall() function will test all the functions in this project. The list
 * of the test functions are extracted from functions below in this file matching
 * the pattern "int test...()". Each test function will test one function for a
 * perticular algorithm, returning PASSED(1) or FAILED(0). At last, a report will
 * be printed after everything is finished.
 */
void testall()
{
    /* the num is defined in "Test.h", which is the number of test functions. */
    int* result = (int*)malloc_s(num * sizeof(int));  /* record the test results */

    /* the tests is defined in "Test.h", which is an array of test functions.*/
    for(int i = 0; i < num; i ++)
    {
        printf("===================================================\n");
        printf("* * * * * * * * * *testing No.%2d* * * * * * * * * *\n", i);
        printf("Function name: %s\n", names[i]);
        result[i] = tests[i]();
        printf("===================================================\n\n");
    }

    for (int i = 0; i < num; i++)
    {
        printf("%s No.%02d %s\n", result[i] ? " [passed]" : "*[failed]", i, names[i]);
    }

    free(result);
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

static double (*integrals[])(double) = {integral1, integral2};

static double (*integrands[])(double) = {integrand1, integrand2};

static int integratecount = 2;

/**
 * @brief test adaptive Simpson method
 * @returns 
 * 
 * 
 */
int testAdaptiveSimpson()
{
    double a = 0;
    double b = 3;
    double TOL = 1e-4;

    double res;
    double ans;
    double (*integrand)(double);
    double (*integral)(double);
    for(int i = 0; i < integratecount; i ++)
    {
        integrand = integrands[i];
        integral = integrals[i];
        res = AdaptiveSimpsonInt(integrand, a, b, TOL);
        ans = integral(b) - integral(a);
        printf("T=%.12f\nI=%.12f\ndelta=%.12f\n", res, ans, fabs(res - ans));
        if(fabs(res - ans) > TOL)
        {
            return FAILED;
        }
    }

    return PASSED;
}

/**
 * @brief 
 * @returns 
 * 
 * 
 */
int testRomberg()
{
    int N = 100;
    double a = 0;
    double b = 1.5;
    double TOL = 1e-10;

    double res;
    double ans;
    double (*integrand)(double);
    double (*integral)(double);
    for(int i = 0; i < integratecount; i ++)
    {
        integrand = integrands[i];
        integral = integrals[i];
        res = RombergInt(integrand, a, b, N, TOL);
        ans = integral(b) - integral(a);
        printf("T=%.12f\nI=%.12f\ndelta=%.12f\n", res, ans, fabs(res - ans));
        if(fabs(res - ans) > TOL)
        {
            return FAILED;
        }
    }

    return PASSED;
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

/**
 * @brief 
 * @returns 
 * 
 * 
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
 * 
 * 
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
 * 
 * 
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
 * 
 * 
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
 * 
 * 
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


/**
 * testing solving ODEs
 */

static double ODEf1(double t, double y)
{
    return (y * y + y) / t;
}

static double ODEy1(double t)
{
    return 2 * t / (1 - 2 * t);
}

static double ODEf2(double t, double y)
{
    return -y + t + 1;
}

static double ODEy2(double t)
{
    return exp(-t) + t;
}

static double ODEf3(double t, double y)
{
    return - y + t * t + 1;
}

static double ODEy3(double t)
{
    return - 2 / exp(t) + 3 - 2 * t + t * t;
}

//static double (*ODEfs[])(double, double) = {ODEf1, ODEf2, ODEf3};
//
//static double (*ODEys[])(double) = {ODEy1, ODEy2, ODEy3};
//
//static int ODEfcount = 3;

/**
 * @brief 
 * @param a 
 * @param b 
 * @param N 
 * @param ODEf 
 * @param ODEy 
 * @param eps 
 * @returns 
 * 
 * 
 */
static int _testClassicRK(double a, double b, int N, double ODEf(double, double), double ODEy(double), double eps)
{
    double t, real_y;
    double y0 = ODEy(a);

    double *result = ClassicRungeKutta(ODEf, a, b, y0, N);
    printf("Testing Classic Runge-Kutta method\n");
    for(int i = 0; i < N + 1; i ++)
    {
        t = a + i * (b - a) / N;
        real_y = ODEy(t);
        printf("%lf\t%lf\n", result[i], real_y);
        if (fabs(result[i] - real_y) > eps)
        {
            fprintf(stderr, "Classic Runge-Kutta method accuracy not enough\n");
            return FAILED;
        }
    }
    return PASSED;
}

/**
 * @brief 
 * @returns 
 * 
 * 
 */
int testClassicRK()
{
    double eps = 0.01;
    return _testClassicRK(1, 3, 4, ODEf1, ODEy1, eps)
        & _testClassicRK(1, 2, 5, ODEf2, ODEy2, eps)
        & _testClassicRK(0, 1, 5, ODEf3, ODEy3, eps);
}


/**
 * @brief 
 * @param a 
 * @param b 
 * @param N 
 * @param ODEf 
 * @param ODEy 
 * @param eps 
 * @returns 
 * 
 * 
 */
static int _testAdamsPECE(double a, double b, int N, double ODEf(double, double), double ODEy(double), double eps)
{
    double y0 = ODEy(a);
    double dy0 = ODEf(a, y0);
    double *ans = (double*)malloc_s((N + 1) * sizeof(double));
    for(int i = 0; i < N + 1; i++)
    {
        ans[i] = ODEy(a + (b - a) * i / N);
    }

    double *result = AdamsPECE(ODEf, a, b, dy0, y0, N);
    printf("Testing Adams PECE\n");
    for(int i = 0; i < N + 1; i ++)
    {
        printf("%lf\t%lf\n", result[i], ans[i]);
        if(fabs(result[i] - ans[i]) > eps)
        {
            fprintf(stderr, "Adams PECE method accuracy not enough\n");
            return FAILED;
        }
    }
    return PASSED;
}

/**
 * @brief 
 * @returns 
 * 
 * 
 */
int testAdamsPECE()
{
    double eps = 1e-3;
    return _testAdamsPECE(0, 1, 10, ODEf2, ODEy2, eps)
        & _testAdamsPECE(1, 2, 20, ODEf1, ODEy1, eps)
        & _testAdamsPECE(0, 1, 10, ODEf3, ODEy3, eps);
}

/**
 * @brief 
 * @param a 
 * @param b 
 * @param TOL 
 * @param hmax 
 * @param hmin 
 * @param ODEf 
 * @param ODEy 
 * @returns 
 * 
 * 
 */
static int _testRKF(double a, double b, double TOL, double hmax, double hmin, double ODEf(double, double), double ODEy(double))
{
    double y0 = ODEy(a);
    double *result;

    printf("Testing RKF Method\n");
    result = RKF78(ODEf, a, b, y0, TOL, hmax, hmin);
    if(result == NULL)
    {
        fprintf(stderr, "RKF method got NULL result.\n");
        return FAILED;
    }
    printf("t\t\th\t\ty\t\treal_y\n");
    for(int i = 0; i < result[0]; i ++)
    {
        for(int j = 0; j < 3; j++)
        {
            printf("%lf\t", result[i * 3 + j + 1]);
        }
        printf("%lf\n", ODEy(result[i * 3 + 1]));

        // Check the results
        if (fabs(result[i * 3 + 3] - ODEy(result[i * 3 + 1])) > 1e-5)
        {
            return FAILED;
        }
    }
    return PASSED;
}

/**
 * @brief 
 * @returns 
 * 
 * 
 */
int testRKF()
{
    return _testRKF(0, 1, 1e-5, 1e-2, 1e-6, ODEf3, ODEy3);
}

static double *SODEf1(double t, double *y)
{
    double *f = (double*)malloc_s(2 * sizeof(double));
    f[0] = -4 * y[0] - 2 * y[1] + cos(t) + 4 * sin(t);
    f[1] = 3 * y[0] + y[1] - 3 * sin(t);
    return f;
}

static double *SODEy1(double x)
{
	double *y = (double*)malloc_s(2 * sizeof(double));
    y[0] = exp(-2 * x) * (-2 + 2 * exp(x) + exp(2 * x) * sin(x));
    y[1] = -exp(-2 * x) * (3 * exp(x) - 2);
    return y;
}

/**
 * @brief 
 * @returns 
 * 
 * 
 */
int testSODERungeKutta()
{
    double *(*f)(double,double*) = SODEf1;
    double *(*y)(double) = SODEy1;
    double a = 0;
    double b = 1;
    double y0[] = {0, -1};
    //int N = 10;
    int m = 2;
    double *t;
    double **res;
    int steps = SODERKF(&t, &res, f, y0, a, b, m, 0.1, 1e-8, 1e4, 1e-4, 13);
    //SODERungeKutta(f, a, b, y0, m, N);

    

    printf("Testing Runge-Kutta Method for a System of ODEs\n");
    printf("%8s%16s%16s%16s%16s\n", "t", "result1", "y1", "result2", "y2");
    for(int i = 0; i < steps; i ++)
    {
		double *ans = y(t[i]);
        printf("%8.2lf", t[i]);
        for(int j = 0 ; j < m; j++)
        {
            printf("%16lf%16lf", res[i][j], ans[j]);
        }
        printf("\n");
    }
    free(t);
    for (int i = 0; i < steps; i++)
	{
		free(res[i]);
	}
	free(res);
	
    return PASSED;
}


/**
 * testing orthogonal decomposition
 */

/**
 * @brief 
 * @returns 
 * 
 * 
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

/**
 * @brief 
 * @returns 
 * 
 * 
 */
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

/**
 * @brief 
 * @returns 
 * 
 * 
 */
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

/**
 * @brief 
 * @returns 
 * 
 * 
 */
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

/**
 * @brief 
 * @returns 
 * 
 * 
 */
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

/**
 * @brief 
 * @returns 
 * 
 * 
 */
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

/**
 * @brief 
 * @returns 
 * 
 * 
 */
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

/**
 * @brief 
 * @returns 
 * 
 * 
 */
int testMuller()
{
    return PASSED;
}
