/**
 * @file Array.c
 * @brief Array manipulation functions
 */

#include "NR.h"
#include "NRprivate.h"

/*--------------------------- one dimension array ----------------------------*/
/**
 * @brief create a one dim array with length n
 * @param N the length of the array
 * @returns a pointer to the array
 */
double *newArray1d(int N)
{
    double *array = (double*)malloc_s(N * sizeof(double));
    return array;
}

/**
 * @brief free the memory of 1d array
 * @param array the 1d array need to be deleted
 */
void delArray1d(double *array)
{
    free(array);
}

/**
 * @brief resize the array
 * @param array array
 * @param M the new size, and it should be larger than the origin length
 * @returns a new array with length M,
 */
double *Array1dResize(double *array, int M)
{
    return (double*)realloc_s((void*)array, M * sizeof(double));
}

/*--------------------------- two dimension array ----------------------------*/
/**
 * @brief create a two dim array with shape N1 * N2
 * @param N1 the length of first axis
 * @param N2 the length of second axis
 * @returns a 2-rank pointer to an array of N1 pointers each of which pointing
 * to an array of N2 double
 */
double **newArray2d(int N1, int N2)
{
    double **array = (double**)malloc_s(N1 * sizeof(double*));
    for(int i = 0; i < N1; i++)
        array[i] = (double*)malloc_s(N2 * sizeof(double));
    return array;
}

/**
 * @brief free the memory of a 2d array
 * @param array the 1d array need to be deleted
 * @param N1 the length of the first axis of array
 */
void delArray2d(double **array, int N1)
{
    for(int i = 0; i < N1; i++)
        free(array[i]);
    free(array);
}

/*------------------------------- utils --------------------------------------*/
/**
 * @brief create the augmented matrix of two matrices with the same row numbers.
 * @param mat1 the first matrix with N rows and M1 columns in 1d form.
 * @param mat2 the second matrix with N rows and M2 columns in 1d form.
 * @param N the row number of two matrices
 * @param M1 the col number of the first matrix
 * @param M2 the col number of the second matrix
 * @returns the augmented matrix of the two matrices in 2d array form
 */
double **AugmentedMatrix(double *mat1, double *mat2, int N, int M1, int M2)
{
    double **AM = newArray2d(N, M1 + M2);
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < M1; j++)
        {
            AM[i][j] = mat1[i * M1 + j];
        }
        for(int j = 0; j < M2; j++)
        {
            AM[i][j + M1] = mat2[i * M2 + j];
        }
    }
    return AM;
}
