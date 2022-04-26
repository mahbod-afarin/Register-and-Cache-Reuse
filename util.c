/**
 * @file util.c
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
#include <sys/time.h>
#include "util.h"

int print_matrix(const double *A, const int m, const int n)
{
	int i;
	printf("[");
	for (i = 0; A + i && i < m * n; i++)
	{

		if ((i + 1) % n == 0)
			printf("%5.2f ", A[i]);
		else
			printf("%5.2f, ", A[i]);
		if ((i + 1) % n == 0)
		{
			if (i + 1 < m * n)
				printf(";\n");
		}
	}
	printf("]\n");

    if (i != m * n) return -1;
    return 0;
}

int randomize_matrix(double *A, const int m, const int n)
{
	srand(time(NULL));
	int i, j;
	for (i = 0; i < m; i++)
	{
		for (j = 0; A + i * n + j && j < n; j++)
		{
			A[i * n + j] = (double)(rand() % 100) + 0.01 * (rand() % 100);
			if ( (rand() % 2) == 0 )
			{
				A[i * n + j] *= -1.;
			}
		}
        if (j != n) return -1;
	}
    return 0;
}

double get_sec()
{
	struct timeval time;
	gettimeofday(&time, NULL);
	return (time.tv_sec + 1e-6 * time.tv_usec);
}


int matrix_copy(double *C, const double *D, const int m, const int n)
{
	int i;

	for (i = 0; C + i && D + i && i < m * n; i++)
	{
		C[i] = D[i];
	}
    if (i != m * n) 
    {
        if (!(C + i)) 
        {
            printf("invalid memory access at %d on input 1, return -1\n", i);
            return -1;
        }
        if (!(D + i))
        {
            printf("invalid memory access at %d on input 2, return -2\n", i);
            return -2;
        }
    }

    return 0;
}

int verify_matrix(const double *C, const double *D, const int m, const int n)
{
	int i;
	double diff = 0.;
    /* we assume input matrices are perfectly initialized */
	for (i = 0; i < m * n; i++)
	{
		diff = fabs(C[i] - D[i]);
        if (diff > 1e-3) break;
	}

	if (diff > 1e-3) 
    {
//        printf("incorrect. bias = %lf\n", diff);
        return -1;
    }
    else
    {
//        printf("correct!\n");
        return 0;
    }
    
}

void *time_measurement(gemm func, char *funcName, double *A, double *B, double *C, int n, double * time) {
    double t0, t1;

	t0 = get_sec();
    func(A, B, C, n);
    t1 = get_sec();
	*time = t1 - t0;

}

void *time_measurement_block(gemm_block func, char *funcName, double *A, double *B, double *C, int n, int b, double * time) {
    double t0, t1;

	t0 = get_sec();
    func(A, B, C, n, b);
    t1 = get_sec();
	*time = t1 - t0;

}