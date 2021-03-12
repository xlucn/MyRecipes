/**
 * @file testSolve.c
 * @brief testSolve.c
 */

#include <stdio.h>
#include <math.h>
#include "NR.h"
#include "Test.h"
#include "constants.h"

static double f1(double x) {return x * x * x - x - 1;}
static double f2(double x) {return (x + 2) * x * x - 4;}
static double g2(double x) {return x - f2(x) / (3 * x + 4) / x;}
static double g3(double x) {return sqrt(4 / (2 + x));}
static double f4(double x) {return ((x + 2) * x + 10) * x - 20;}
static double df4(double x) {return (3 * x + 4) * x + 10;}
static double f5(double x) {return (2 * x * x - 5) * x - 1;}

/**
 * @brief Solve test unit
 */
typedef struct SolveTest{
    double (*f)(double); /**< empty */
    double (*g)(double); /**< empty */
    double (*df)(double); /**< empty */
    double root; /**< empty */
    double x1; /**< empty */
    double x2; /**< empty */
    double x3; /**< empty */
    double TOL; /**< empty */
}SolveTest;

static SolveTest solvetests[] = {
/*   f   g     df    root              x1   x2   x3   TOL   */
    {f1, NULL, NULL, 1.32471795724475, 1.0, 2.0, NAN, 1e-13},
    {f2, g2,   NULL, 1.130395435,      1.0, NAN, NAN, 1e-8 },
    {f2, g3,   NULL, 1.130395435,      1.5, NAN, NAN, 1e-8 },
    {f4, NULL, df4,  1.368808108,      1.0, NAN, NAN, 1e-8 },
    {f5, NULL, NULL, 1.672981647854,   2.0, 1.0, 1.5, 1e-11},
    {NULL}
};

static int _testSolve(double (*f)(SolveTest), SolveTest t)
{
    double p = f(t);
    printf("%lf\n", p);
    if(fabs(p - t.root) > t.TOL)
    {
        return FAILED;
    }
    return PASSED;
}

/*--------------------------- Bisection Method -------------------------------*/
static double _testBisection(SolveTest t)
{
    return Bisection(t.f, t.x1, t.x2, t.TOL);
}

/**
 * @brief testBisection
 * @return integer if test passed
 */
int testBisection()
{
    for(int i = 0; solvetests[i].f; i++) if(!isnan(solvetests[i].x2))
    {
        if(_testSolve(_testBisection, solvetests[i]) == FAILED)
            return FAILED;
    }
    return PASSED;
}

/*---------------------------- Picard Recurtion-------------------------------*/
static double _testPicardRecurtion(SolveTest t)
{
    return PicardIteration(t.g, t.x1, t.TOL);
}

/**
 * @brief testPicardRecurtion
 * @return integer if test passed
 */
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
    return SteffensenIteration(t.g, t.x1, t.TOL);
}

/**
 * @brief Steffensen
 * @return integer if test passed
 */
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
    return NewtonMethod(t.f, t.df, t.x1, t.TOL);
}

/**
 * @brief NewtonMethod
 * @return integer if test passed
 */
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
    return SecentMethod(t.f, t.x1, t.x2, t.TOL);
}

/**
 * @brief Secent
 * @return integer if test passed
 */
int testSecent()
{
    for(int i = 0; solvetests[i].f; i++) if(!isnan(solvetests[i].x2))
    {
        if(_testSolve(_testSecent, solvetests[i]) == FAILED)
            return FAILED;
    }
    return PASSED;
}

/*---------------------------- Muller Method ---------------------------------*/
static double _testMuller(SolveTest t)
{
    return MullerMethod(t.f, t.x1, t.x2, t.x3, t.TOL);
}

/**
 * @brief Muller
 * @return integer if test passed
 */
int testMuller()
{
    for(int i = 0; solvetests[i].f; i++) if(!isnan(solvetests[i].x3))
    {
        if(_testSolve(_testMuller, solvetests[i]) == FAILED)
            return FAILED;
    }
    return PASSED;
}
