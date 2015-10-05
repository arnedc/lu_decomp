#include <stdio.h>
#include <mkl.h>
#include <mkl_cblas.h>
#include <stdlib.h>
#include <stddef.h>
#include "matrix.h"

dense_matrix multiply_blas_matrix(int n, dense_matrix matA, dense_matrix matB)
{
  double *matprod=(dense_matrix) malloc(n*n*sizeof(double));
  if (matprod==NULL){
    printf("Unable to create dense matrix (function: multiply_blas_matrix)");
    return NULL;
  }
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,n,n,n,1,matA,n,matB,n,0,matprod,n);
  return matprod;
}