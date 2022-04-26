/**
 * @file reg_reuse.c
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
    printf("\n*********** Register Reuse ***********\n");
    int i, j;
    
    char *func_name_reg[] = {"dgemm0", "dgemm1", "dgemm2", "dgemm3"};
    int matrix_dimension[] = {66, 126, 258, 510, 1026, 2046};
    int m = sizeof(matrix_dimension) / sizeof(matrix_dimension[0]);

    int func_nums = sizeof(func_name_reg) / sizeof(func_name_reg[0]);
    double t0, t1;
    double *result = (double *)malloc(sizeof(double) * m * func_nums * 2);    
    gemm algorithm_ptr[] = {&dgemm0, &dgemm1, &dgemm2, &dgemm3};

    for (j = 0; j < m; j++)
    {
        int curr_dim = matrix_dimension[j];
        double *A = (double *)malloc(sizeof(double) * curr_dim * curr_dim);
        double *B = (double *)malloc(sizeof(double) * curr_dim * curr_dim);
        if ( randomize_matrix(A, curr_dim, curr_dim) ) return -1;
        if ( randomize_matrix(B, curr_dim, curr_dim) ) return -1;
        double *C[func_nums];
        double *C_verify[func_nums];
        for (i = 0; i < func_nums; i++)
        {
            C[i] = (double *)malloc(sizeof(double) * curr_dim * curr_dim);
            C_verify[i] = (double *)malloc(sizeof(double) * curr_dim * curr_dim);
            if ( randomize_matrix(C[i], curr_dim, curr_dim) ) return -1;
            if ( matrix_copy(C_verify[i], C[i], curr_dim, curr_dim) ) return -1;
        }

        for (i = 0; i < func_nums; i++)
        {
            time_measurement(algorithm_ptr[i], func_name_reg[i], A, B, C[i], curr_dim, &result[j * func_nums * 2 + i * 2]);
            t0 = get_sec();
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, curr_dim, curr_dim, curr_dim, 1., A, curr_dim, B, curr_dim, 1., C_verify[i], curr_dim);
            t1 = get_sec();
            result[j * func_nums * 2 + 2 * i + 1] = t1 - t0;
            if ( verify_matrix(C[i], C_verify[i], curr_dim, curr_dim) )
            {
                printf("error detected at function %s, matrix size = %d.\n", func_name_reg[i], curr_dim);
            }
        }

        free(A);
        free(B);
        for (i = 0; i < func_nums; i++)
        {
            free(C[i]);
            free(C_verify[i]);
        }
    }

    for (i = 0; i < func_nums; i++)
    {
        for (j = 0; j < m; j++)
        {
            printf("%s - %4d: elapsed time is %8.5f second(s).\n", func_name_reg[i], matrix_dimension[j], result[2 * j * func_nums + 2 * i]);
        }
    }

    free(result);
    return 0;
}