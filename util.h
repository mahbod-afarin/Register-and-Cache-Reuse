#ifndef _UTIL_H_
#define _UTIL_H_

#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef void (*gemm)(double *, double *, double *, int);
typedef void (*gemm_block)(double *, double *, double *, int, ...);

int print_matrix(const double *A, const int m, const int n);
int randomize_matrix(double *A, const int m, const int n);
double get_sec();
int matrix_copy(double *C, const double *D, const int m, const int n);
int verify_matrix(const double *C, const double *D, const int m, const int n);
void *time_measurement(gemm func, char *funcName, double *A, double *B, double *C, int n, double * time);
void *time_measurement_block(gemm_block func, char *funcName, double *A, double *B, double *C, int n, int b, double * time);

#endif // _UTIL_H_