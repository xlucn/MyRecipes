//LuXu
//Solving an equation
#include <stdlib.h>
#include <math.h>
#include "NumericalRecipes.h"
#include "LibFunction.h"

double Bisection(double f(double), double a, double b, double eps)
{
    double p;
    while(fabs(b - a) > eps)
    {
        p = (a + b) / 2;
        (f(a) * f(p) < 0) ? (b = p) : (a = p);
    }
    return p;
}

double PicardIteration(double g(double), double x, double eps)
{
    double p = g(x);
    while(fabs(p - x) > eps)
    {
        x = p;
        p = g(x);
    }
    return p;
}

double SteffensenIteration(double g(double), double x, double eps)
{
    double y = g(x);
    double z = g(y);
    double p = x - (y - x) * (y - x) / (z - 2 * y + x);
    while(fabs(p - x) > eps)
    {
        x = p;
        y = g(x);
        z = g(y);
        p = x - (y - x) * (y - x) / (z - 2 * y + x);
    }
    return p;
}

double NewtonMethod(double f(double), double df(double), double x, double eps)
{
    double p = x - f(x) / df(x);
    while(fabs(p - x) > eps)
    {
        x = p;
        p = x - f(x) / df(x);
    }
    return p;
}

double SecentMethod(double f(double), double x0, double x1, double eps)
{
    double y0 = f(x0);
    double y1 = f(x1);
    double p = x1 - y1 * (x1 - x0) / (y1 - y0);
    while(fabs(p - x1) > eps)
    {
        x0 = x1;
        x1 = p;
        y0 = y1;
        y1= f(p);
        p = x1 - y1 * (x1 - x0) / (y1 - y0);
    }
    return p;
}

double MullerMethod(double f(double), double x0, double x1, double x2, double eps)
{
    double b;
    double d;
    double e;
    double h;

    double f0 = f(x0);
    double f1 = f(x1);
    double f2 = f(x2);
    double f01 = (f1 - f0) / (x1 - x0);
    double f12 = (f2 - f1) / (x2 - x1);
    double f012  = (f12 - f01) / (x2 - x0);

    for(int i = 0; i < 1000; i++)
    {
        b = f12 + f012 * (x2 - x1);
        d = sqrt(b * b - 4 * f2 * f012);
        e = fabs(b - d) > fabs(b + d) ? b + d : b - d;
        h = -2 * f2 / e;
        if (fabs(h) < eps)
        {
            break;
        }
        x0 = x1;
        x1 = x2;
        x2 += h;
        f0 = f1;
        f1 = f2;
        f2 = f(x2);
        f01 = (f1 - f0) / (x1 - x0);
        f12 = (f2 - f1) / (x2 - x1);
        f012  = (f12 - f01) / (x2 - x0);
    }

    return x2;
}
