#include "NR.h"
#include "NRprivate.h"

double *newArray1d(int N)
{
    double *array = (double*)malloc_s(N * sizeof(double));
    return array;
}

void delArray1d(double *array)
{
    free(array);
}

double **newArray2d(int N1, int N2)
{
    double **array = (double**)malloc_s(N1 * sizeof(double*));
    for(int i = 0; i < N1; i++)
        array[i] = (double*)malloc_s(N2 * sizeof(double));
    return array;
}

/**
 * @brief 
 * @param array 
 * @param N1 
 */
void delArray2d(double **array, int N1)
{
    for(int i = 0; i < N1; i++)
        free(array[i]);
    free(array);
}

/**
 * @brief create the augmented matrix of two matrices with the same row numbers.
 * @param mat1 the first matrix with N rows and M1 columns in 1d form.
 * @param mat2 the second matrix with N rows and M2 columns in 1d form.
 * @param N the row number of two matrices
 * @param M1 the col number of the first matrix
 * @param M2 the col number of the second matrix
 * @returns the augmented matrix of the two matrices
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
