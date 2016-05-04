/**
 * @file NR.h
 * @note encoding : utf-8
 * @author Lu Xu(oliver_lew@outlook.com)
 * @brief 这里是根据林成森《数值计算方法》一书编写的计算程序。
 */
#ifndef _NR_H_
#define _NR_H_


/*****************************************************************************
 *********************LinearEquations*****************************************
 *****************************************************************************/
/**
 * @brief Chasing method for solving tridiagonal equations.
 * @param N the order of the matrix.
 * @param d the main diagonal,
 * @param c the line above d,
 * @param a the line below d,
 * @param b the constant vector.
 * @return the solution of the equation *x or NULL if the eqation has no solution.
 */
double *Chasing(int N, double *d, double *c, double *a, double *b);
/**
 * @brief Gaussian elimination method to solve linear equation in the form of Ax=b.
 * @param N The rank of the matrix
 * @param A The coefficient matrix
 * @param b The constant vector
 */
double *GaussEli(int N, double **A, double *b);
/**
 * @brief Gauss elimination with partial pivoting
 * @param N the rank of the equation
 * @param a the augmented matrix
 */
double *GaussEliPP(int N, double **a);
/**
 * @brief Gauss elimination with partial pivoting proportionally
 * @param N the rank of the equation
 * @param a the augmented matrix
 */
double *GaussEliPPP(int N, double **a);
/**
 * @brief Gauss Jordan elimination method to solve a system of linear equations
 * @param N the rank of the equation
 * @param a the augmented matrix
 */
double *GaussJordanEli(int N, double **a);



/*****************************************************************************
 *********************Integral**********************************************
 *****************************************************************************/
/**
 * @brief Trapezoidal integration
 */
double TrapezoidalInt(double(*f)(double), double a, double b);
/**
 * @brief Simpson integration
 */
double SimpsonInt(double(*f)(double), double a, double b);
/**
 * @brief Composite Simpson integration，N is the number of the subintervals
 */
double CompositeSimpsonInt(double(*f)(double), double a, double b, int N);
/**
 * @brief Romberg integration，
 * @param N 所设置的最多序列数,
 * @param [a,b] the interval,
 * @param eps the precision,
 * @param f the integration function
 */
double RombergInt(double(*f)(double), double a, double b, int N, double eps);
/**
 * @brief Adaptive Simpson intergral, using recursion.
 */
double AdaptiveSimpsonInt(double(*f)(double), double a, double b, double TOL);



/*****************************************************************************
 *********************ODE*****************************************************
 *****************************************************************************/
/**
 * @brief a struct type to contain the solution of a system of ODE integration
 * 
 * This struct type contains three menbers: steps, record the steps used in the
 * integration. It is useful in variable step size methods. t, the pointer to the
 * independent virable list. y, the pointer to the dependent variables list.
 */
typedef struct _SODEsol{
	int step; /**< the steps used in the integration */
	double *t; /**< the pointer to the independent virable list */
	double **y; /**< the pointer to the dependent variables list */
}SODEsol;

/**
 * @brief Euler method to solve initial value problem(IVP) of ODE
 * @param y0 initial value,
 * @param f derivative function. dy/dt=f(t,y).
 * @param a lower limit of interval
 * @param b upper limit of interval
 * @param N number of subintervals
 */
double* Euler(double(*f)(double, double), double a, double b, double y0, int N);
/**
 * @brief Improved Euler Method
 */
double* ImprovedEuler(double(*f)(double, double), double a, double b, double y0, int N);
/**
 * @brief Midpoint method.
 */
double* MID(double(*f)(double, double), double a, double b, double y0, int N);
/**
 * @brief Heun Method
 */
double* Heun(double(*f)(double, double), double a, double b, double y0, int N);
/**
 * @brief Three stage Heun method.
 */
double *ThreeStageHeunMethod(double(*f)(double, double), double a, double b, double y0, int N);
/**
 * @brief Three Stage Runge-Kutta Method
 */
double *ThreeStageRungeKuttaMathod(double(*f)(double, double), double a, double b, double y0, int N);
/**
 * @brief Classic Runge-Kutta Method
 */
double *ClassicRungeKutta(double(*f)(double, double), double a, double b, double y0, int N);
/**
 * @brief Runge-Kutta-Fehlberg Method
 */
double *RKF78(double(*f)(double,double), double a, double b, double y0, double TOL, double hmax, double hmin);
double *RKF45(double(*f)(double,double), double a, double b, double y0, double TOL, double hmax, double hmin);
/**
 * @brief One-step Adams correlation PECE method. Use classic Runge-Kutta method for the initial value.
 * @param f right function of ODE
 * @param a lower limit of interval
 * @param b upper limit of interval
 * @param dy0 initial value of derivative of f
 * @param y0 initial value of f
 * @param N number of subintervals
 * @return array of ys
 */
double *AdamsPECE(double(*f)(double, double), double a, double b, double dy0, double y0, int N);
/**
 * @brief Classic Runge-Kutta Method to solve a System of ODEs.
 * @param m the number of ODEs,
 * @param N the numebr of steps,
 * @param (a, b) the interval,
 * @param y0 the initial values,
 * @param f point to an array of functions.
 * @return a 2D array of values of all functions in all steps.
 */
double **SODERungeKutta(double (*f[])(double, double*), double a, double b, double *y0, int m, int N);
/**
 * @brief RKF method to solve a system of ODEs
 * @param t the address of a pointer to double type, this function will modify
 *  the pointer to point to a one-dimension array.
 * @param y the address of a 2-rank pointer to double type, this function will
 *  modify the pointer to point to a two-dimension array.
 * @param f an array of pointers to right-functions
 * @param y0 an array of initial values
 * @param a,b interval
 * @param m number of equations
 * @param h0 the initial step size
 * @param TOL the required tolerance
 * @param hmin, hmax the minimum and maximum step size allowed
 * @param n temporarily use this number to indicate which method to use
 * @return the number of steps use by RKF method to solve the equations, the
 *  pointer *t will be a array of time in every step, and the pointer *y will be
 *  a 2-rank array of all function values in all steps.
 * Example Usage:
 * @code
 *
 * @endcode
 */
SODEsol SODERKF(double *(*f)(double, double*), double *y0,
     double a, double b, int m, double h0, double TOL, double hmax, double hmin, int n);


/*****************************************************************************
 *********************Interpolation*****************************************
 *****************************************************************************/
/**
 * @brief Hermite polynomial Interpolation
 * @param N number of Interpolation points,
 * @param a[N] the Interpolation points,
 * @param x variable
 * @param f[N] function,
 * @param df[N] derivative of f with respect to variable x,
 * @return approximate value using Hermite polynomial
 */
double Hermite(int N, double *a, double *f, double *df, double x);
/**
 * @brief 自然三次样条插值函数
 * @param N 分段区间数，即插值点数量为N + 1；
 * @param a 插值点数组指针；
 * @param f 函数；
 * @param x 变量值；
 * @param ddf_a second order derivative at a
 * @param ddf_b second order derivative at b
 * @return 样条插值函数值
 */
double NatureCubicSplineIpl(double(*f)(double), double ddf_a, double ddf_b, double x, int N, double *a);
/**
 * @brief 完备三次样条插值函数
 * @param N 分段区间数，即插值点数量为N + 1；
 * @param a 插值点数组指针；
 * @param f 函数；
 * @param x 变量值；
 * @param df_a derivative at a
 * @param df_b derivative at b
 * @return 样条插值函数值
 */
double CompleteCubicSplineIpl(double(*f)(double), double df_a, double df_b, double x, int N, double *a);
/**
 * @brief Lagrange三次样条插值函数
 * @param N 分段区间数，即插值点数量为N + 1；
 * @param a 插值点数组指针；
 * @param f 函数；
 * @param x 变量值；
 * @return 返回样条插值函数值
 */
double LagrangeCubicSplineIpl(double(*f)(double), double x, int N, double *a);



/*****************************************************************************
 *********************Basic*************************************************
 *****************************************************************************/
/**
 * @brief Divided difference of function f on nodes x[n]
 * f[x] = f(x),
 * f[x1, x2] = (f(x2)-f(x1))/(x2-x1),
 * ......
 * f[x1, x2, ... , xn] = (f[x2, x3, ..., xn]-f[x1, x2, ..., x(n-1)])/(xn-x1)
 */
double DividedDiff(double(*f)(double), double *x, int N);
/**
 * @brief All the Divided differences of a array up to order k
 */
double **FullDividedDiff(double(*f)(double), double *x, int N, int k);
/**
 * @brief return a matrix of the divided differnces. d[i][j] is f[xi, ..., xj]
 */
double **DividedDiffMatrix(double(*f)(double), double *x, int N);
/**
 * @brief Lagrange polynomial
 * 输入N和插值点向量a[N]，以及变量值x，返回拉格朗日基本多项式的值
 */
double *LagrangePoly(double *a, double x, int N);
/**
 * @brief Chebyshev polynomial(first kind)
 */
double Chebyshev(int n, double x);



/******************************************************************************
**********************LeastSq************************************************
******************************************************************************/
/**
 * @brief solve the leastsquare solution for a system of linear equations.
 */
double *LeastSquare(int m, int n, double **A, double *b);
/**
 * @brief Gram-Schmidt method of POD(proper orthogonal decomposition).
 * Rank of The m*n(m>n) matrix need to be n.
 * return a three demension array which contains two 2-d arrays [Q,R].
 */
double ***GramSchmidtQR(int m, int n, double **A);
/**
 * @brief Improved version of Gram-Schmidt method.
 */
double ***ImprovedGramSchmidtQR(int m, int n, double **A);



/*****************************************************************************
 *********************Solve*************************************************
 *****************************************************************************/
/**
 * @brief Bisection method to find the root of a equation
 */
double Bisection(double(*f)(double), double a, double b, double eps);
/**
 * @brief Use integration to solve equations in the form of x = g(x)
 */
double PicardIteration(double g(double), double x, double eps);
/**
 * @brief Accelerate iteration of PicardIteration
 */
double SteffensenIteration(double g(double), double x, double eps);
/**
 * @brief Newton method (or Newton-Raphson Method) to find the root of a equation.
 * df is the derivative function of f(x)
 */
double NewtonMethod(double(*f)(double), double df(double), double x, double eps);
/**
 * @brief Secent method to solve a equation, x0 and x1 is two initial points you need.
 */
double SecentMethod(double(*f)(double), double x0, double x1, double eps);
/**
 * @brief Muller method
 */
double MullerMethod(double(*f)(double), double x0, double x1, double x2, double eps);

#endif // _NR_H_
