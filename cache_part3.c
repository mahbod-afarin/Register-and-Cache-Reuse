/**
 * @file cache_part3.c
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
    printf("\n*********** Cache reuse ***********\n");
    int i, j;
    int b;
    char *func_name_cache[] = {"ijk ", "jik ", "kij ", "ikj ", "jki ", "kji ",
                               "bijk", "bjik", "bkij", "bikj", "bjki", "bkji"};
    int matrix_dim = 2048;
    int block_size=64;
    int method_nums = sizeof(func_name_cache) / sizeof(func_name_cache[0]);
    double t0, t1;
    double *result = (double *)malloc(sizeof(double) * method_nums * 2);
    gemm_block algorithm_ptr[] = {&ijk, &jik, &kij, &ikj, &jki, &kji, &bijk, &bjik, &bkij, &bikj, &bjki, &bkji};

    double *A = (double *)malloc(sizeof(double) * matrix_dim * matrix_dim);
    double *B = (double *)malloc(sizeof(double) * matrix_dim * matrix_dim);
    if ( randomize_matrix(A, matrix_dim, matrix_dim) ) return -1;
    if ( randomize_matrix(B, matrix_dim, matrix_dim) ) return -1;
    double *C[method_nums];
    double *C_verify[method_nums];
    for (i = 0; i < method_nums; i++)
    {
        C[i] = (double *)malloc(sizeof(double) * matrix_dim * matrix_dim);
        C_verify[i] = (double *)malloc(sizeof(double) * matrix_dim * matrix_dim);
        if ( randomize_matrix(C[i], matrix_dim, matrix_dim) ) return -1;
        if ( matrix_copy(C_verify[i], C[i], matrix_dim, matrix_dim) ) return -1;
    }

    for (i = 0; i < method_nums; i++)
    {
        time_measurement_block(algorithm_ptr[i], func_name_cache[i], A, B, C[i], matrix_dim, block_size, &result[i * 2]);
        t0 = get_sec();
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, matrix_dim, matrix_dim, matrix_dim, 1., A, matrix_dim, B, matrix_dim, 1., C_verify[i], matrix_dim);
        t1 = get_sec();
        result[2 * i + 1] = t1 - t0;
        if (verify_matrix(C[i], C_verify[i], matrix_dim, matrix_dim))
        {
            printf("error detected at function %s, matrix size = %d.\n", func_name_cache[i], matrix_dim);
        }
    }

    free(A);
    free(B);
    for (i = 0; i < method_nums; i++)
    {
        free(C[i]);
        free(C_verify[i]);
    }

    for (i = 0; i < method_nums; i++)
    {
        printf("%s: elapsed time is %8.5f second(s).\n", func_name_cache[i], result[2 * i]);
    }

    free(result);
    return 0;
}
