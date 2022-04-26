/**
 * @file cache_part4.c
 * author Yujia Zhai (yzhai015@ucr.edu)
 * @brief 
 * @version 0.1
 * @date 2019-10-08
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mkl.h"

#include "util.h"
#include "mygemm.h"

int main(int argc, char *argv[])
{
    int i, j;
    int b;

    int matrix_dim = 2046;
    int block_size= 66;
    double t0, t1, elapsed_time, MKL_elapsed_time;
    double *A = (double *)malloc(sizeof(double) * matrix_dim * matrix_dim);
    double *B = (double *)malloc(sizeof(double) * matrix_dim * matrix_dim);
    double *C = (double *)malloc(sizeof(double) * matrix_dim * matrix_dim);
    double *C_verify = (double *)malloc(sizeof(double) * matrix_dim * matrix_dim);

    if ( randomize_matrix(A, matrix_dim, matrix_dim) ) return -1;
    if ( randomize_matrix(B, matrix_dim, matrix_dim) ) return -1;
    if ( randomize_matrix(C, matrix_dim, matrix_dim) ) return -1;
    if ( matrix_copy(C_verify, C, matrix_dim, matrix_dim) ) return -1;

    t0 = get_sec();
    optimal(A, B, C, matrix_dim,block_size);
    t1 = get_sec();
    elapsed_time = t1 - t0;

    t0 = get_sec();
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, matrix_dim, matrix_dim, matrix_dim, 1., A, matrix_dim, B, matrix_dim, 1., C_verify, matrix_dim);
    t1 = get_sec();
    MKL_elapsed_time = t1 - t0;
    if (verify_matrix(C, C_verify, matrix_dim, matrix_dim))
    {
        printf("error detected at part4, matrix size = %d.\n", matrix_dim);
    }

    free(A);
    free(B);
    free(C);
    free(C_verify);

    printf("elapsed time is %8.5f second(s).\n", elapsed_time);

    return 0;
}
