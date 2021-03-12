/**
 * @file Test.c
 * @brief test functions
 * @note if a function is declared as "int test...()", it will be collected by
 *       GenerateTest.py and regarded as a test function.
 */

#include <stdio.h>
#include "Test.h"

/**
 * @brief Test all the functions in this project
 *
 * This will test all the functions in this project. The test
 * functions list is extracted from other source files in this folder. All the
 * function names matching the pattern "int test...()" will be included. Each
 * test function will test one function for a perticular algorithm, returning
 * PASSED(0) or FAILED(1). At last, a report will be printed after everything
 * is finished.
 *
 * @return if all tests passed
 */
int main()
{
    /*
     * Macros defined in "Test.h":
     * "FUNC_COUNT" : the number of test functions.
     * "FUNC_ARRAY" : an array of test functions.
     * "NAME_ARRAY" : an array of test functions' names.
     */
    int (*tests[FUNC_COUNT])() = FUNC_ARRAY;
    char *names[FUNC_COUNT] = NAME_ARRAY;
    int result[FUNC_COUNT], result_all = PASSED;

    for(int i = 0; i < FUNC_COUNT; i++)
    {
        printf("\n===================================================\n");
        printf("* * * * * * * * * *testing No.%2d* * * * * * * * * *\n", i);
        printf("Function name: %s\n", names[i]);
        result[i] = tests[i]();
        printf("\n===================================================\n\n");
    }

    /* summary */
    for (int i = 0; i < FUNC_COUNT; i++)
    {
        printf("%s No.%02d %s\n", (result[i] == PASSED) ?
            "\033[36m[passed]\033[0m" : "\033[31m[failed]\033[0m", i, names[i]);
        if (result[i] == FAILED)
            result_all = FAILED;
    }

    return result_all;
}
