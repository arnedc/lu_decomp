#include <stdlib.h>
#include "matrix.h"
#include <stdio.h>
#include <mkl_cblas.h>
#include <mkl_lapack.h>
#include <mkl_lapacke.h>
#include <math.h>
#include <string.h>

int lapack_symm_solve(int n,dense_matrix matA, int *P, double *matB)
{
    return LAPACKE_dpptrs(LAPACK_ROW_MAJOR,'L',n,1,matA,matB,1);
}

int lapack_solve(int n, dense_matrix matA, int *P, double *matB)
{
    return LAPACKE_dgetrs(LAPACK_ROW_MAJOR,'L',n,1,matA,n,P,matB,1);
}

double *dense_solve(int n, dense_matrix matA, int *P, double *matB)
{
    int i,j,rn;
    double temp, *sol;
    sol=calloc(n,sizeof(double));
    for (i=0;i<n;++i) {
        if (*(P+i)==-1) {
            *(sol+i)=NAN;
            continue;
        }
        if (*(P+i)!=i) {
            rn=*(P+i);
            temp=*(matB+i);
            *(matB+i)=*(matB + rn);
            *(matB + rn)=temp;
        }
        *(sol+i) = *(matB+i);
        for (j=i-1;j>=0;--j) {
            if (!(isnan(*(sol+j))))
                *(sol+i) -= *(sol+j) * *(matA+i*n+j);
        }
    }
    for (i=n-1;i>=0;--i) {
        if (!(isnan(*(sol+i)))) {
            for (j=i+1;j<n;++j) {
                if (!(isnan(*(sol+j))))
                    *(sol+i) -= *(sol+j) * *(matA+i*n+j);
            }
            *(sol+i) /= *(matA+i*n+i);
        }
    }
    return sol;
}

double *dense_cholsolve(int n, dense_matrix matA, double *matB)
{
  int i,j,rn;
    double temp, *sol;
    sol=calloc(n,sizeof(double));
    for (i=0;i<n;++i) {
        if (isnan(*(matA+(i*i+3*i)/2))) {
            *(sol+i)=NAN;
            continue;
        }
        *(sol+i) = *(matB+i);
        for (j=i-1;j>=0;--j) {
            if (!(isnan(*(sol+j))))
                *(sol+i) -= *(sol+j) * *(matA+(i*i+i)/2+j);
        }
        *(sol+i) /= *(matA+(i*i+3*i)/2);
    }
    for (i=n-1;i>=0;--i) {
        if (!(isnan(*(sol+i)))) {
            for (j=i+1;j<n;++j) {
                if (!(isnan(*(sol+j))))
                    *(sol+i) -= *(sol+j) * *(matA+(j*j+j)/2+i);
            }
            *(sol+i) /= *(matA+(i*i+3*i)/2);
        }
    }
    return sol;
}

void dense_trisolve(int n, dense_matrix matA, double matB)
{
  
}

void dense_RFP_solve(int n, dense_matrix matA, double *matB)
{
  if(n%2==0){
    cblas_dtrsm(CblasRowMajor,CblasLeft,CblasLower,CblasNoTrans,CblasNonUnit,n/2,1,1,matA+n/2,n/2,matB,1);
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,n/2,1,n/2,-1,matA+(n/2+1)*n/2,n/2,matB,1,1,matB+n/2,1);
    cblas_dtrsm(CblasRowMajor,CblasLeft,CblasUpper,CblasTrans,CblasNonUnit,n/2,1,1,matA,n/2,matB+n/2,1);
    cblas_dtrsm(CblasRowMajor,CblasLeft,CblasUpper,CblasNoTrans,CblasNonUnit,n/2,1,1,matA,n/2,matB+n/2,1);
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,n/2,1,n/2,-1,matA+(n/2+1)*n/2,n/2,matB+n/2,1,1,matB,1);
    cblas_dtrsm(CblasRowMajor,CblasLeft,CblasLower,CblasTrans,CblasNonUnit,n/2,1,1,matA+n/2,n/2,matB,1);
  }
  else{
    cblas_dtrsm(CblasRowMajor,CblasLeft,CblasLower,CblasNoTrans,CblasNonUnit,(n+1)/2,1,1,matA,(n+1)/2,matB,1);
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,(n-1)/2,1,(n+1)/2,-1,matA+(n+1)/2*(n+1)/2,(n+1)/2,matB,1,1,matB+(n+1)/2,1);
    cblas_dtrsm(CblasRowMajor,CblasLeft,CblasUpper,CblasTrans,CblasNonUnit,(n-1)/2,1,1,matA+1,(n+1)/2,matB+(n+1)/2,1);
/*    cblas_dtrsm(CblasRowMajor,CblasLeft,CblasUpper,CblasNoTrans,CblasNonUnit,(n-1)/2,1,1,matA+1,(n+1)/2,matB+(n+1)/2,1);
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,(n+1)/2,1,(n-1)/2,-1,matA+(n+1)/2*(n+1)/2,(n+1)/2,matB+(n+1)/2,1,1,matB,1);
    cblas_dtrsm(CblasRowMajor,CblasLeft,CblasLower,CblasTrans,CblasNonUnit,(n+1)/2,1,1,matA,(n+1)/2,matB,1);
  */}
}