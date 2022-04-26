#include "mygemm.h"

/**
 * 
 * Implement all functions here in this file.
 * Do NOT change input parameters and return type.
 * 
 **/

//Register Reuse part 1
void dgemm0(const double* A, const double* B, double* C, const int n)
{
    int i = 0, j = 0, k = 0;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            for (k = 0; k < n; k++)
            {
                C[i * n + j] += A[i * n + k] * B[k * n + j];
            }
        }
    }
}

void dgemm1(const double* A, const double* B, double* C, const int n)
{
    int i = 0, j = 0, k = 0;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            register double C_Temp = C[i * n + j];
            for (k = 0; k < n; k++)
            {
                C_Temp += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = C_Temp;
        }
    }
}
//Register Reuse part 1 End

//Register Reuse part 2
void dgemm2(const double* A, const double* B, double* C, const int n)
{
    int i = 0, j = 0, k = 0;
    for (i = 0; i < n; i += 2)
    {
        for (j = 0; j < n; j += 2)
        {
            register double C1 = C[i * n + j];
            register double C2 = i < (n - 1) ? C[(i + 1) * n + j] : 0;
            register double C3 = j < (n - 1) ? C[i * n + (j + 1)] : 0;
            register double C4 = (i < (n - 1)) && (j < (n - 1))? C[(i + 1) * n + (j + 1)] : 0;
            for (k = 0; k < n; k += 2)
            {
                register double A1 = A[i * n + k];
                register double A2 = i < (n - 1) ? A[(i + 1) * n + k] : 0;
                register double A3 = k < (n - 1) ? A[i * n + (k + 1)] : 0;
                register double A4 = (i < (n - 1)) && (k < (n - 1)) ? A[(i + 1) * n + (k + 1)] : 0;
                register double B1 = B[k * n + j];
                register double B2 = k < (n - 1) ? B[(k + 1) * n + j] : 0;
                register double B3 = j < (n - 1) ? B[k * n + (j + 1)] : 0;
                register double B4 = (k < (n - 1)) && (j < (n - 1)) ? B[(k + 1) * n + (j + 1)] : 0;

                C1 = C1 + A1 * B1 + A3 * B2;
                C2 = C2 + A2 * B1 + A4 * B2;
                C3 = C3 + A1 * B3 + A3 * B4;
                C4 = C4 + A2 * B3 + A4 * B4;
            }

            C[i * n + j] = C1;
            C[(i + 1) * n + j] = C2;
            C[i * n + (j + 1)] = C3;
            C[(i + 1) * n + (j + 1)] = C4;
        }
    }
}

//Register Reuse part 2 End

//Register Reuse part 3
void dgemm3(const double* A, const double* B, double* C, const int n)
{
    int i = 0, j = 0, k = 0;
    for (i = 0; i < n; i += 3)
    {
        for (j = 0; j < n; j += 4)
        {
            register double C1 = C[i * n + j];
            register double C2 = i < (n - 1) ? C[(i + 1) * n + j] : 0;
            register double C3 = i < (n - 2) ? C[(i + 2) * n + j] : 0;
            register double C4 = j < (n - 1) ? C[i * n + (j + 1)] : 0;
            register double C5 = i < (n - 1) && j < (n - 1) ? C[(i + 1) * n + (j + 1)] : 0;
            register double C6 = i < (n - 2) && j < (n - 1) ? C[(i + 2) * n + (j + 1)] : 0;
            register double C7 = j < (n - 2) ? C[i * n + (j + 2)] : 0;
            register double C8 = i < (n - 1) && j < (n - 2) ? C[(i + 1) * n + (j + 2)] : 0;
            register double C9 = i < (n - 2) && j < (n - 2) ? C[(i + 2) * n + (j + 2)] : 0;
            register double C10 = j < (n - 3) ? C[i * n + (j + 3)] : 0;
            register double C11 = i < (n - 1) && j < (n - 3) ? C[(i + 1) * n + (j + 3)] : 0;
            register double C12 = i < (n - 2) && j < (n - 3) ? C[(i + 2) * n + (j + 3)] : 0;
            for (k = 0; k < n; k++)
            {
                register double A1 = A[i * n + k];
                register double A2 = i < (n - 1) ? A[(i + 1) * n + k] : 0;
                register double A3 = i < (n - 2) ? A[(i + 2) * n + k] : 0;
                register double B1 = B[k * n + j];

                C1 = C1 + A1 * B1;
                C2 = C2 + A2 * B1;
                C3 = C3 + A3 * B1;
                B1 = j < (n - 1) ? B[k * n + (j + 1)] : 0;
                C4 = C4 + A1 * B1;
                C5 = C5 + A2 * B1;
                C6 = C6 + A3 * B1;
                B1 = j < (n - 2) ? B[k * n + (j + 2)] : 0;
                C7 = C7 + A1 * B1;
                C8 = C8 + A2 * B1;
                C9 = C9 + A3 * B1;
                B1 = j < (n - 3) ? B[k * n + (j + 3)] : 0;
                C10 = C10 +  A1 * B1;
                C11 = C11 + A2 * B1;
                C12 = C12 + A3 * B1;
            }

            C[i * n + j] = C1;
            if (i < (n - 1)) C[(i + 1) * n + j] = C2;
            if (i < (n - 2)) C[(i + 2) * n + j] = C3;
            if (j < (n - 1)) C[i * n + (j + 1)] = C4;
            if (i < (n - 1) && j < (n - 1)) C[(i + 1) * n + (j + 1)] = C5;
            if (i < (n - 2) && j < (n - 1)) C[(i + 2) * n + (j + 1)] = C6;
            if (j < (n - 2)) C[i * n + (j + 2)] = C7;
            if (i < (n - 1) && j < (n - 2)) C[(i + 1) * n + (j + 2)] = C8;
            if (i < (n - 2) && j < (n - 2)) C[(i + 2) * n + (j + 2)] = C9;
            if (j < (n - 3)) C[i * n + (j + 3)] = C10;
            if (i < (n - 1) && j < (n - 3)) C[(i + 1) * n + (j + 3)] = C11;
            if (i < (n - 2) && j < (n - 3)) C[(i + 2) * n + (j + 3)] = C12;
        }
    }
}
//Register Reuse part 3 End

//Cache Reuse part 3
void ijk(const double *A, const double *B, double *C, const int n) 
{
    int i = 0, j = 0, k = 0;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            register double Cij = C[i * n + j];
            for (k = 0; k < n; k++)
            {
                Cij += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = Cij;
        }
    }
}

void bijk(const double *A, const double *B, double *C, const int n, const int b) 
{
    int i = 0, j = 0, k = 0;
    int kBLK = 0, jBLK = 0, iBLK = 0;
    for(i=0; i<n; i+=b)
    {
		for(j=0; j<n; j+=b)	
        {
			for(k=0; k<n; k+=b)
			{
				for(iBLK = i; iBLK < i + b && iBLK < n; iBLK++)
                {
					for(jBLK = j; jBLK < j + b && jBLK < n; jBLK++)
					{
						register double temp = C[iBLK * n + jBLK];
						for(kBLK = k; kBLK < k + b && kBLK < n; kBLK++)
                        {
							temp += A[iBLK * n + kBLK] * B[kBLK * n + jBLK];
                        }
						C[iBLK * n + jBLK] = temp;
					}
                }
			}
        }
    }
}

void jik(const double *A, const double *B, double *C, const int n) 
{
    int i = 0, j = 0, k = 0;
    for (j = 0; j < n; j++)
    {
        for (i = 0; i < n; i++)
        {
            register double Cij = C[i * n + j];
            for (k = 0; k < n; k++)
            {
                Cij += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = Cij;
        }
    }
}

void bjik(const double *A, const double *B, double *C, const int n, const int b) 
{
    int i = 0, j = 0, k = 0;
    int kBLK = 0, jBLK = 0, iBLK = 0;
    for( j=0;j<n;j+=b)
    {
		for(i=0;i<n;i+=b)	
        {
			for(k=0;k<n;k+=b)
			{
				for(jBLK=j;jBLK<j+b && jBLK<n;jBLK++)
                {
					for(iBLK=i;iBLK<i+b && iBLK<n;iBLK++)
					{
						register double temp=C[iBLK*n+jBLK];
						for(kBLK=k;kBLK<k+b && kBLK<n;kBLK++)
                        {
							temp+=A[iBLK*n+kBLK]*B[kBLK*n+jBLK];
                        }
						C[iBLK*n+jBLK]=temp;
					}
                }    
			}
        }
    }    
}

void kij(const double *A, const double *B, double *C, const int n) 
{
    int i = 0, j = 0, k = 0;
    for (k = 0; k < n; k++)
    {
        for (i = 0; i < n; i++)
        {
            register double Aik = A[i * n + k];
            for (j = 0; j < n; j++)
            {
                C[i * n + j] += Aik * B[k * n + j];
            }
        }
    }
}

void bkij(const double *A, const double *B, double *C, const int n, const int b) 
{
    int i = 0, j = 0, k = 0;
    int kBLK = 0, jBLK = 0, iBLK = 0;
    for(k=0;k<n;k+=b)
    {
		for(i=0;i<n;i+=b)	
        {
			for(j=0;j<n;j+=b)
			{
				for(kBLK=k;kBLK<k+b && kBLK<n;kBLK++)
                {
					for(iBLK=i;iBLK<i+b && iBLK<n;iBLK++)
					{
						register double temp=A[iBLK*n+kBLK];
						for(jBLK=j;jBLK<j+b && jBLK<n;jBLK++)
                        {
							C[iBLK*n+jBLK]+=temp*B[kBLK*n+jBLK];
                        }    
					}
                }    
			}
        }
    }        
}


void ikj(const double *A, const double *B, double *C, const int n) 
{
    int i = 0, j = 0, k = 0;
    for (i = 0; i < n; i++)
    {
        for (k = 0; k < n; k++)
        {
            register double Aik = A[i * n + k];
            for (j = 0; j < n; j++)
            {
                C[i * n + j] += Aik * B[k * n + j];
            }
        }
    }
}

void bikj(const double *A, const double *B, double *C, const int n, const int b) 
{
    int i = 0, j = 0, k = 0;
    int kBLK = 0, jBLK = 0, iBLK = 0;
    for(i=0;i<n;i+=b)
    {
		for(k=0;k<n;k+=b)
        {	
			for(j=0;j<n;j+=b)
			{
				for(iBLK=i;iBLK<i+b && iBLK<n;iBLK++)
                {
					for(kBLK=k;kBLK<k+b && kBLK<n;kBLK++)
					{
						register double temp=A[iBLK*n+kBLK];
						for(jBLK=j;jBLK<j+b && jBLK<n;jBLK++)
                        {
							C[iBLK*n+jBLK]+=temp*B[kBLK*n+jBLK];
                        }    
					}
                }    
			}
        }
    }        
}

void jki(const double *A, const double *B, double *C, const int n) 
{
    int i = 0, j = 0, k = 0;
    for (j = 0; j < n; j++)
    {
        for (k = 0; k < n; k++)
        {
            register double Bkj = B[k * n + j];
            for (i = 0; i < n; i++)
            {
                C[i * n + j] += A[i * n + k] * Bkj;
            }
        }
    }
}

void bjki(const double *A, const double *B, double *C, const int n, const int b) 
{
    int i = 0, j = 0, k = 0;
    int kBLK = 0, jBLK = 0, iBLK = 0;
    for(j=0;j<n;j+=b)
    {
		for(k=0;k<n;k+=b)	
        {
			for(i=0;i<n;i+=b)
			{
				for(jBLK=j;jBLK<j+b && jBLK<n;jBLK++)
                {
					for(kBLK=k;kBLK<k+b && kBLK<n;kBLK++)
					{
						register double temp=B[kBLK*n+jBLK];
						for( iBLK=i;iBLK<i+b && iBLK<n;iBLK++)
                        {
							C[iBLK*n+jBLK]+=temp*A[iBLK*n+kBLK];
                        }    
					}
                }    
			}
        }
    }        
}

void kji(const double *A, const double *B, double *C, const int n) 
{
    int i = 0, j = 0, k = 0;
    for (k = 0; k < n; k++)
    {
        for (j = 0; j < n; j++)
        {
            register double Bkj = B[k * n + j];
            for (i = 0; i < n; i++)
            {
                C[i * n + j] += A[i * n + k] * Bkj;
            }
        }
    }
}

void bkji(const double *A, const double *B, double *C, const int n, const int b) 
{
    int i = 0, j = 0, k = 0;
    int kBLK = 0, jBLK = 0, iBLK = 0;
    for(k=0;k<n;k+=b)
    {
		for(j=0;j<n;j+=b)	
        {
			for(i=0;i<n;i+=b)
			{
				for(kBLK=k;kBLK<k+b && kBLK<n;kBLK++)
                {
					for(jBLK=j;jBLK<j+b && jBLK<n;jBLK++)
					{
						register double temp=B[kBLK*n+jBLK];
						for(iBLK=i;iBLK<i+b && iBLK<n;iBLK++)
                        {
							C[iBLK*n+jBLK]+=temp*A[iBLK*n+kBLK];
                        }    
					}
                }    
			}
        }
    }        
}
//Cache Reuse part 3 End 

//Cache Reuse part 4
void optimal(const double* A, const double* B, double *C, const int n, const int b)
{
    int i = 0, j = 0, k = 0;
    int iBLK = 0, jBLK = 0, kBLK = 0;
    for (i = 0; i < n; i += b)
    {
        for (j = 0; j < n; j += b)
        {
            for (k = 0; k < n; k += b)
            {
                for (iBLK = i; iBLK < (i + b > n? n : (i + b)); iBLK += 3)
                {
                    for (jBLK = j; jBLK < (j + b > n? n : (j + b)); jBLK += 3)
                    {
                        register double C1 = C[iBLK * n + jBLK];
                        register double C2 = C[(iBLK + 1) * n + jBLK];
                        register double C3 = C[(iBLK + 2) * n + jBLK];
                        register double C4 = C[iBLK * n + (jBLK + 1)];
                        register double C5 = C[(iBLK + 1) * n + (jBLK + 1)];
                        register double C6 = C[(iBLK + 2) * n + (jBLK + 1)];
                        register double C7 = C[iBLK * n + (jBLK + 2)];
                        register double C8 = C[(iBLK + 1) * n + (jBLK + 2)];
                        register double C9 = C[(iBLK + 2) * n + (jBLK + 2)];

                        for (kBLK = k; kBLK < (k + b > n? n : (k + b)); kBLK++)
                        {
                            register double B0 =  B[kBLK * n + jBLK];
                            register double A0 = A[iBLK * n + kBLK];
                            register double A1 = A[(iBLK + 1) * n + kBLK];
                            register double A2 = A[(iBLK + 2) * n + kBLK];
                            C1 = C1 + A0 * B0;
                            C2 = C2 + A1 * B0;
                            C3 = C3 + A2 * B0;
                            B0 = B[kBLK * n + (jBLK + 1)];
                            C4 = C4 + A0 * B0;
                            C5 = C5 + A1 * B0;
                            C6 = C6 + A2 * B0;
                            B0 = B[kBLK * n + (jBLK + 2)];
                            C7 = C7 + A0 * B0;
                            C8 = C8 + A1 * B0;
                            C9 = C9 + A2 * B0;

                        }
                        C[iBLK * n + jBLK] = C1;
                        C[(iBLK + 1) * n + jBLK] = C2;
                        C[(iBLK + 2) * n + jBLK] = C3;
                        C[iBLK * n + (jBLK + 1)] = C4;
                        C[(iBLK + 1) * n + (jBLK + 1)] = C5;
                        C[(iBLK + 2) * n + (jBLK + 1)] = C6;
                        C[iBLK * n + (jBLK + 2)] = C7;
                        C[(iBLK + 1) * n + (jBLK + 2)] = C8;
                        C[(iBLK + 2) * n + (jBLK + 2)] = C9;
                    }
                }
            }
        }
    }
}