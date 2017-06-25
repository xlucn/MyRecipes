/** @file testSolve.c */
#include <stdio.h>
#include <math.h>
#include "NR.h"
#include "Test.h"
#include "constants.h"

/** tolerance allowed for the functions */
#define TOL 1e-8

static double f1(double x) {return x * x * x - x - 1;}
static double f2(double x) {return (x + 2) * x * x - 4;}
static double g2(double x) {return x - f2(x) / (3 * x + 4) / x;}
static double g3(double x) {return sqrt(4 / (2 + x));}
static double f4(double x) {return ((x + 2) * x + 10) * x - 20;}
static double df4(double x) {return (3 * x + 4) * x + 10;}
static double f5(double x) {return (2 * x * x - 5) * x - 1;}

typedef struct SolveTest{
    double (*f)(double);
    double (*g)(double);
    double (*df)(double);
    double root;
    double x1;
    double x2;
    double x3;
}SolveTest;

static SolveTest solvetests[] = {
    {f1,   NULL, NULL, 1.32471795724475, 1.0, 2.0      , FLOAT_NAN},
    {f2,   g2,   NULL, 1.130395435,      1.0, FLOAT_NAN, FLOAT_NAN},
    {f2,   g3,   NULL, 1.130395435,      1.5, FLOAT_NAN, FLOAT_NAN},
    {f4,   NULL, df4,  1.368808108,      1.0, FLOAT_NAN, FLOAT_NAN},
    {f5,   NULL, NULL, 1.672981647854,   2.0, 1.0      , 1.5      },
    {NULL}
};

static int _testSolve(double (*f)(SolveTest), SolveTest t)
{
    double p = f(t);
    printf("%lf\n", p);
    if(fabs(p - t.root) > TOL)
    {
        return FAILED;
    }
    return PASSED;
}

/*--------------------------- Bisection Method -------------------------------*/
static double _testBisection(SolveTest t)
{
    return Bisection(t.f, t.x1, t.x2, TOL);
}

int testBisection()
{
    for(int i = 0; solvetests[i].f; i++) if(solvetests[i].x2 != FLOAT_NAN)
    {
        if(_testSolve(_testBisection, solvetests[i]) == FAILED)
            return FAILED;
    }
    return PASSED;
}

/*---------------------------- Picard Recurtion-------------------------------*/
static double _testPicardRecurtion(SolveTest t)
{
    return PicardIteration(t.g, t.x1, TOL);
}

int testPicardRecurtion()
{
    for(int i = 0; solvetests[i].f; i++) if(solvetests[i].g != NULL)
    {
        if(_testSolve(_testPicardRecurtion, solvetests[i]) == FAILED)
            return FAILED;
    }
    return PASSED;
}

/*------------------------ Steffensen iteration-------------------------------*/
static double _testSteffensen(SolveTest t)
{
    return SteffensenIteration(t.g, t.x1, TOL);
}

int testSteffensen()
{
    for(int i = 0; solvetests[i].f; i++) if(solvetests[i].g != NULL)
    {
        if(_testSolve(_testSteffensen, solvetests[i]) == FAILED)
            return FAILED;
    }
    return PASSED;
}

/*---------------------------- Newton Method ---------------------------------*/
static double _testNewtonMethod(SolveTest t)
{
    return NewtonMethod(f4, df4, t.x1, TOL);
}

int testNewtonMethod()
{
    for(int i = 0; solvetests[i].f; i++) if(solvetests[i].df != NULL)
    {
        if(_testSolve(_testNewtonMethod, solvetests[i]) == FAILED)
            return FAILED;
    }
    return PASSED;
}

/*---------------------------- Secent Method ---------------------------------*/
static double _testSecent(SolveTest t)
{
    return SecentMethod(t.f, t.x1, t.x2, TOL);
}

int testSecent()
{
    for(int i = 0; solvetests[i].f; i++) if(solvetests[i].x2 != FLOAT_NAN)
    {
        if(_testSolve(_testSecent, solvetests[i]) == FAILED)
            return FAILED;
    }
    return PASSED;
}

/*---------------------------- Muller Method ---------------------------------*/
static double _testMuller(SolveTest t)
{
    return MullerMethod(t.f, t.x1, t.x2, t.x3, TOL);
}

int testMuller()
{
    for(int i = 0; solvetests[i].f; i++) if(solvetests[i].x3 != FLOAT_NAN)
    {
        if(_testSolve(_testMuller, solvetests[i]) == FAILED)
            return FAILED;
    }
    return PASSED;
}
