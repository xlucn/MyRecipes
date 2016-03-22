//encoding : utf-8
//Lu Xu
//这里是根据林成森《数值计算方法》一书编写的计算程序。
#ifndef _NR_H_
#define _NR_H_


/*****************************************************************************
 *********************LinearEquations*****************************************
 *****************************************************************************/
/**
 * Chasing method for solving tridiagonal equations.
 * @param N the order of the matrix.
 * @param d the main diagonal,
 * @param c the line above d,
 * @param a the line below d,
 * @param d the constant vector.
 * @return the solution of the equation *x or NULL if the eqation has no solution.
 */
double *Chasing(int N, double *d, double *c, double *a, double *b);
/**
 * Gaussian elimination method to solve linear equation in the form of Ax=b.
 */
double *GaussEli(int N, double **A, double *b);
/**
 * Gauss elimination with partial pivoting
 */
double *GaussEliPP(int N, double **a);
/**
 * Gauss elimination with partial pivoting proportionally
 */
double *GaussEliPPP(int N, double **a);
/**
 * Gauss Jordan elimination method to solve a system of linear equations
 */
double *GaussJordanEli(int N, double **a);



/*****************************************************************************
 *********************Integral.c**********************************************
 *****************************************************************************/
/**
 * Trapezoidal integration
 */
double TrapezoidalInt(double f(double), double a, double b);
/**
 * Simpson integration
 */
double SimpsonInt(double f(double), double a, double b);
/**
 * Composite Simpson integration，N is the number of the subintervals
 */
double CompositeSimpsonInt(double f(double), double a, double b, int N);
/**
 * Romberg integration，
 * @param N 所设置的最多序列数,
 * @param [a,b] the interval,
 * @param eps the precision,
 * @param f the integration function
 */
double RombergInt(double f(double), double a, double b, int N, double eps);
/**
 * Adaptive Simpson intergral, using recursion.
 */
double AdaptiveSimpsonInt(double f(double), double a, double b, double TOL);



/*****************************************************************************
 *********************ODE.c***************************************************
 *****************************************************************************/
/**
 * Euler method to solve initial value problem(IVP) of ODE
 * @param y0 initial value,
 * @param f derivative function. dy/dt=f(t,y).
 */
double* Euler(double f(double, double), double a, double b, double y0, int N);
/**
 * Improved Euler Method
 */
double* ImprovedEuler(double f(double, double), double a, double b, double y0, int N);
/**
 * Midpoint method.
 */
double* MID(double f(double, double), double a, double b, double y0, int N);
/**
 * Heun Method
 */
double* Heun(double f(double, double), double a, double b, double y0, int N);
/**
 * Three stage Heun method.
 */
double *ThreeStageHeunMethod(double f(double, double), double a, double b, double y0, int N);
/**
 * Three Stage Runge-Kutta Method
 */
double *ThreeStageRungeKuttaMathod(double f(double, double), double a, double b, double y0, int N);
/**
 * Classic Runge-Kutta Method
 */
double *ClassicRungeKutta(double f(double, double), double a, double b, double y0, int N);
/**
 * Runge-Kutta-Fehlberg Method
 */
double *RKF78(double f(double,double), double a, double b, double y0, double TOL, double hmax, double hmin);
double *RKF45(double f(double,double), double a, double b, double y0, double TOL, double hmax, double hmin);
/**
 * Adams显式和隐式方法的PECE模式校正方法，这里k=1，用经典Runge-Kutta方法提供初值
 */
double *AdamsPECE(double f(double, double), double a, double b, double dy0, double y0, int N);
/**
 * Classic Runge-Kutta Method to solve a System of ODEs.
 * @param m the number of ODEs,
 * @param N the numebr of steps,
 * @param (a, b) the interval,
 * @param y0 the initial values,
 * @param f point to an array of functions.
 * @return a 2D array of values of all functions in all steps.
 */
double **SODERungeKutta(double (*f[])(double, double*), double a, double b, double *y0, int m, int N);
/**
 * RKF method to solve system of ODEs
 */
 int SODERKF(double** t, double*** y, double (**f)(double, double*), double *y0,
     double a, double b, int m, double h0, double TOL, double hmax, double hmin, int n);


/*****************************************************************************
 *********************Interpolation.c*****************************************
 *****************************************************************************/
/**
 * Hermite polynomial Interpolation
 * @param N number of Interpolation points,
 * @param a[N] the Interpolation points,
 * @param f[N] function,
 * @param df[N] derivative of f with respect to variable x,
 * @return approximate value using Hermite polynomial
 */
double Hermite(int N, double *a, double *f, double *df, double x);
/**
 * 自然三次样条插值函数
 * @param N 分段区间数，即插值点数量为N + 1；
 * @param a 插值点数组指针；
 * @param f 函数；
 * @param x 变量值；
 * @param ddf_a ddf_b second order derivative at a and b respectively;
 * @return 样条插值函数值
 */
double NatureCubicSplineIpl(double f(double), double ddf_a, double ddf_b, double x, int N, double *a);
/**
 * 完备三次样条插值函数
 * @param N 分段区间数，即插值点数量为N + 1；
 * @param a 插值点数组指针；
 * @param f 函数；
 * @param x 变量值；
 * @param df_a, df_b两端导数值；
 * @return 样条插值函数值
 */
double CompleteCubicSplineIpl(double f(double), double df_a, double df_b, double x, int N, double *a);
/**
* Lagrange三次样条插值函数
* @param N 分段区间数，即插值点数量为N + 1；
* @param a 插值点数组指针；
* @param f 函数；
* @param x 变量值；
* @return 返回样条插值函数值
*/
double LagrangeCubicSplineIpl(double f(double), double x, int N, double *a);



/*****************************************************************************
 *********************Basic.c*************************************************
 *****************************************************************************/
/**
 * Divided difference of function f on nodes x[n]
 * f[x] = f(x),
 * f[x1, x2] = (f(x2)-f(x1))/(x2-x1),
 * ......
 * f[x1, x2, ... , xn] = (f[x2, x3, ..., xn]-f[x1, x2, ..., x(n-1)])/(xn-x1)
 */
double DividedDiff(double f(double), double *x, int N);
/**
* All the Divided differences of a array up to order k
*/
double **FullDividedDiff(double f(double), double *x, int N, int k);
/**
 * return a matrix of the divided differnces. d[i][j] is f[xi, ..., xj]
 */
double **DividedDiffMatrix(double f(double), double *x, int N);
/**
 * Lagrange polynomial
 * 输入N和插值点向量a[N]，以及变量值x，返回拉格朗日基本多项式的值
 */
double *LagrangePoly(double *a, double x, int N);
/**
 * Chebyshev polynomial(first kind)
 */
double Chebyshev(int n, double x);



/******************************************************************************
**********************LeastSq.c************************************************
******************************************************************************/
/**
 * solve the leastsquare solution for a system of linear equations.
 */
double *LeastSquare(int m, int n, double **A, double *b);
/**
 * Gram-Schmidt method of POD(proper orthogonal decomposition).
 * Rank of The m*n(m>n) matrix need to be n.
 * return a three demension array which contains two 2-d arrays [Q,R].
 */
double ***GramSchmidtQR(int m, int n, double **A);
/**
 * Improved version of Gram-Schmidt method.
 */
double ***ImprovedGramSchmidtQR(int m, int n, double **A);



/*****************************************************************************
 *********************Solve.c*************************************************
 *****************************************************************************/
/**
 * Bisection method to find the root of a equation
 */
double Bisection(double f(double), double a, double b, double eps);
/**
 * Use integration to solve equations in the form of x = g(x)
 */
double PicardIteration(double g(double), double x, double eps);
/**
 * Accelerate iteration of PicardIteration
 */
double SteffensenIteration(double g(double), double x, double eps);
/**
 * Newton method (or Newton-Raphson Method) to find the root of a equation.
 * df is the derivative function of f(x)
 */
double NewtonMethod(double f(double), double df(double), double x, double eps);
/**
 * Secent method to solve a equation, x0 and x1 is two initial points you need.
 */
double SecentMethod(double f(double), double x0, double x1, double eps);
/**
 * Muller method
 */
double MullerMethod(double f(double), double x0, double x1, double x2, double eps);

#endif // _NR_H_
