/**
 * This file is automatically generated by GenerateTest.py, collecting function
 * information from source file "Test.c".
 *
 * If you want to change the content of this file, go to GenerateTest.py and
 * modify the script. Or you can maintain this file by yourself after you delete
 * the related lines in the makefile.
 */

#ifndef _TESTDECLARATION_H_
#define _TESTDECLARATION_H_

#define PASSED 1
#define FAILED 0

extern int (*tests[])();    //the list of test functions which are listed below
extern char *names[];       //the list of test functions' names
extern int num;	            //the number of test functions

//declarations of testfunctions
int testAdamsPECE();
int testAdaptiveSimpson();
int testBisection();
int testChasing();
int testClassicRK();
int testDividedDiff();
int testGaussianEli();
int testGaussianEliPP();
int testHermiteIpl();
int testLeastSq();
int testMuller();
int testNewtonMethod();
int testPicardRecurtion();
int testQR();
int testRKF();
int testRomberg();
int testSODERungeKutta();
int testSecent();
int testSplineIpl();
int testSteffensen();


void testall();
#endif
