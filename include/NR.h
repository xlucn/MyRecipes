/**
 * @file NR.h
 * @note encoding : utf-8
 * @author Lu Xu(oliver_lew@outlook.com)
 * @brief 这里是根据林成森《数值计算方法》一书编写的计算程序。
 */
#ifndef _NR_H_
#define _NR_H_

#include "constants.h"
#include <stdlib.h>

#define malloc
#define realloc
void* malloc_s(size_t size);
void* realloc_s(void* ptr, size_t size);

/*----------------------- Solve linear equations -----------------------------*/
double *Chasing(int N, double *d, double *c, double *a, double *b);
double *GaussEli(int N, double *a, double *b);
double *GaussEliPP(int N, double *a, double *b);
double *GaussEliPPP(int N, double *a, double *b);
double *GaussJordanEli(int N, double **a);

/*--------------------------- Integration ------------------------------------*/
double TrapezoidalInt(double(*f)(double), double a, double b);
double SimpsonInt(double(*f)(double), double a, double b);
double CompositeSimpsonInt(double(*f)(double), double a, double b, int N);
double RombergInt(double(*f)(double), double a, double b, int N, double eps);
double AdaptiveSimpsonInt(double(*f)(double), double a, double b, double TOL);

/*------------------- ODE (ordinary differential equations) ------------------*/
typedef struct _SODEsol *SODEsol;
typedef struct _ODEsol *ODEsol;
struct _SODEsol;
struct _ODEsol;

/* manipulating the structures */
SODEsol newSODEsol(int step, double *t, double **y);
ODEsol newODEsol(int step, double *t, double *y);
int      ODEsolGetStep(ODEsol sol);
double*  ODEsolGetT(ODEsol sol);
double*  ODEsolGetY(ODEsol sol);
void     DisposeODEsol(ODEsol sol);
int      SODEsolGetStep(SODEsol sol);
double*  SODEsolGetT(SODEsol sol);
double** SODEsolGetY(SODEsol sol);
void     DisposeSODEsol(SODEsol sol);

double* Euler(double(*f)(double, double), double a, double b, double y0, int N);
double* ImprovedEuler(double(*f)(double, double), double a, double b, double y0, int N);
double* MID(double(*f)(double, double), double a, double b, double y0, int N);
double* Heun(double(*f)(double, double), double a, double b, double y0, int N);
double* ThreeStageHeunMethod(double(*f)(double, double), double a, double b, double y0, int N);
double* ThreeStageRungeKuttaMathod(double(*f)(double, double), double a, double b, double y0, int N);
double* ClassicRungeKutta(double(*f)(double, double), double a, double b, double y0, int N);
double* RKF78(double(*f)(double,double), double a, double b, double y0, double TOL, double hmax, double hmin);
double* RKF45(double(*f)(double,double), double a, double b, double y0, double TOL, double hmax, double hmin);
double* AdamsPECE(double(*f)(double, double), double a, double b, double dy0, double y0, int N);
double** SODERungeKutta(double (*f[])(double, double*), double a, double b, double *y0, int m, int N);
SODEsol SODERKF(double *(*f)(double, double*), double *y0, double a, double b, int m, double h0, double TOL, double hmax, double hmin, int n);


/*----------------------------- Interpolation --------------------------------*/
double Hermite(int N, double *a, double *f, double *df, double x);
double NatureCubicSplineIpl(double(*f)(double), double ddf_a, double ddf_b, double x, int N, double *a);
double CompleteCubicSplineIpl(double(*f)(double), double df_a, double df_b, double x, int N, double *a);
double LagrangeCubicSplineIpl(double(*f)(double), double x, int N, double *a);

/*------------------------- Basic functions ----------------------------------*/
double DividedDiff(double(*f)(double), double *x, int N);
double **FullDividedDiff(double(*f)(double), double *x, int N, int k);
double **DividedDiffMatrix(double(*f)(double), double *x, int N);
double *LagrangePoly(double *a, double x, int N);
double Chebyshev(int n, double x);

/*------------------------------ Least Square --------------------------------*/
double *LeastSquare(int m, int n, double **A, double *b);
double ***GramSchmidtQR(int m, int n, double **A);
double ***ImprovedGramSchmidtQR(int m, int n, double **A);

/*---------------------------- Solve equations -------------------------------*/
double Bisection(double(*f)(double), double a, double b, double eps);
double PicardIteration(double (*g)(double), double x, double eps);
double SteffensenIteration(double (*g)(double), double x, double eps);
double NewtonMethod(double(*f)(double), double (*df)(double), double x, double eps);
double SecentMethod(double(*f)(double), double x0, double x1, double eps);
double MullerMethod(double(*f)(double), double x0, double x1, double x2, double eps);

#endif // _NR_H_
