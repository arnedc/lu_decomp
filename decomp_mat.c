#include <stdlib.h>
#include "matrix.h"
#include <stdio.h>
#include <mkl_cblas.h>
#include <mkl_lapack.h>
#include <mkl_lapacke.h>
#include <math.h>
#include <string.h>

void dense_LU(int n,dense_matrix matA, dense_tri L, dense_tri U, int *P)
{
    int i,j,m;
    double temp;
    for (i=0;i<n;++i) {
        *(P+i)=i;
        for (j=0;j<n;++j) {
            *(U+j+i*n)=*(matA+j+i*n);
        }
    }
    for (i=0;i<n;++i) {
        *(L+i+i*n)=1;
        for (m=i+1;m<n;++m) {
            if (fabs(*(U+i+m*n)) > fabs(*(U+i+ *(P+i) * n)))
                *(P+i)=m;
        }
        if (*(P+i)!=i) {
            for (j=0;j<n;++j) {
                if (j<i) {
                    temp = *(L+i*n+j);
                    *(L+i*n+j)=*(L + *(P+i) * n + j);
                    *(L + *(P+i) * n + j)=temp;
                }
                else {
                    temp = *(U+i*n+j);
                    *(U+i*n+j)=*(U + *(P+i) * n + j);
                    *(U + *(P+i) * n + j)=temp;
                }
            }
        }
        for (j=i+1;j<n;++j) {
            *(L+i+j*n) = *(U+i+j*n) / *(U+i*n+i);
            for (m=i;m<n;++m) {
                *(U+j*n+m) = *(U+j*n+m) - (*(L+i+j*n) * *(U+i*n+m));
            }
        }
    }
}

void dense_cholesky(int n,dense_matrix matA)
{
    int i,j,k;
    for (i=0;i<n;++i) {
        if (*(matA+(i*i+3*i)/2) > 1e-5) {
            *(matA+(i*i+3*i)/2)=sqrt(*(matA+(i*i+3*i)/2));
            for (j=i+1;j<n;++j) {
                *(matA+(j*j+j)/2+i) = *(matA+(j*j+j)/2+i) / *(matA+(i*i+3*i)/2);
                for (k=i+1;k<j+1;++k) {
                    *(matA+(j*j+j)/2+k) -= *(matA+(j*j+j)/2+i) * *(matA+(k*k+k)/2+i);
                }
            }
        }
        else 
	  *(matA+(i*i+3*i)/2)=NAN;
    }
}

void dense_RFP_cholesky(int n,dense_matrix matA)
{
  int i,*length, *lda;
  length=malloc(sizeof(int));
  lda=malloc(sizeof(int));
  if(n % 2 == 0){
    *lda=n/2;
    for (i=0;i<n/2;++i) {
        if (*(matA+(i+1)*n/2+i) > 1e-5) {
            *(matA+(i+1)*n/2+i)=sqrt(*(matA+(i+1)*n/2+i));
	    *length=n-i-1;
            drscl(length,matA+(i+1)*n/2+i,matA+(i+2)*n/2+i,lda);
	    cblas_dsyr(CblasRowMajor,CblasLower,n/2-i-1,-1,matA+(i+2)*n/2+i,n/2,matA+(i+2)*n/2+i+1,n/2);
	    cblas_dger(CblasRowMajor,n/2,n/2-i-1,-1,matA+(n/2+1)*n/2+i,n/2,matA+(i+2)*n/2+i,n/2,matA+n/2*(n/2+1)+i+1,n/2);
	    cblas_dsyr(CblasRowMajor,CblasUpper,n/2,-1,matA+(n/2+1)*n/2+i,n/2,matA,n/2);
            }
        else 
	  *(matA+(i+1)*n/2+i)=NAN;
    }
 /*   cblas_dtrsm(CblasRowMajor,CblasRight,CblasLower,CblasTrans,CblasNonUnit,n/2,n/2,1,matA+n/2,n/2,matA+(n/2+1)*n/2,n/2);
    cblas_dsyrk(CblasRowMajor,CblasUpper,CblasNoTrans,n/2,n/2,-1,matA+(n/2+1)*n/2,n/2,1,matA,n/2); */   
    *lda=1;
    for (i=0;i<n/2;++i) {
        if (*(matA+i*n/2+i) > 1e-5) {
            *(matA+i*n/2+i)=sqrt(*(matA+i*n/2+i));
	    *length=n/2-i-1;
            drscl(length,matA+i*n/2+i,matA+i*n/2+i+1,lda);
	    cblas_dsyr(CblasRowMajor,CblasUpper,n/2-i-1,-1,matA+i*n/2+i+1,1,matA+(i+1)*n/2+i+1,n/2);
            }
        else 
	  *(matA+i*n/2+i)=NAN;
    }
  }
  else{
    *lda=(n+1)/2;
    for (i=0;i<(n+1)/2;++i) {
        if (*(matA+i*(n+1)/2+i) > 1e-5) {
            *(matA+i*(n+1)/2+i)=sqrt(*(matA+i*(n+1)/2+i));
	    *length=(n+1)/2-i-1;
            drscl(length,matA+i*(n+1)/2+i,matA+(i+1)*(n+1)/2+i,lda);
	    cblas_dsyr(CblasRowMajor,CblasLower,(n+1)/2-i-1,-1,matA+(i+1)*(n+1)/2+i,(n+1)/2,matA+(i+1)*(n+1)/2+i+1,(n+1)/2);
            }
        else 
	  *(matA+i*(n+1)/2+i)=NAN;
    }
    cblas_dtrsm(CblasRowMajor,CblasRight,CblasLower,CblasTrans,CblasNonUnit,(n-1)/2,(n+1)/2,1,matA,(n+1)/2,matA+(n+1)/2*(n+1)/2,(n+1)/2);
    cblas_dsyrk(CblasRowMajor,CblasUpper,CblasNoTrans,(n-1)/2,(n+1)/2,-1,matA+(n+1)/2*(n+1)/2,(n+1)/2,1,matA+1,(n+1)/2);
  *lda=1;
    for (i=0;i<(n-1)/2;++i) {
        if (*(matA+i*(n+1)/2+i+1) > 1e-5) {
            *(matA+i*(n+1)/2+i+1)=sqrt(*(matA+i*(n+1)/2+i+1));
	    *length=(n-1)/2-i;
            drscl(length,matA+i*(n+1)/2+i+1,matA+i*(n+1)/2+i+2,lda);
	    cblas_dsyr(CblasRowMajor,CblasUpper,(n-1)/2-i,-1,matA+i*(n+1)/2+i+2,1,matA+(i+1)*(n+1)/2+i+2,(n+1)/2);
            }
        else 
	  *(matA+i*(n+1)/2+i+1)=NAN;
    }
  }
}

int lapack_LU(int n,dense_matrix matA, int *P)
{
    return LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n,n,matA,n,P);
}

int lapack_symm_LU(int n,dense_matrix matA,int *P)
{
  dense_matrix ARF;
  int info;
  ARF=calloc((n*n+n)/2,sizeof(double));
  if(ARF==NULL){
    printf("Unable to create RFP matrix (function: lapack_symm_lu)");
    return;
  }
  LAPACKE_dtpttf(LAPACK_COL_MAJOR,'N','U',n,matA,ARF);
  info = LAPACKE_dpftrf(LAPACK_COL_MAJOR,'N','U',n,ARF);
  LAPACKE_dtfttp(LAPACK_COL_MAJOR,'N','U',n,ARF,matA);
  return info;

}

void dense_LU_small(int n,int m,int lda, dense_matrix matA, int *P)
{
    int i,j, *length;
    double temp;
    length=malloc(sizeof(int));
    for (i=0;i<n;++i) {
        *(P+i)=i;
    }
    for (i=0;i<n;++i) {
        for (j=i+1;j<m;++j) {
            if (fabs(*(matA+i+j*lda)) > fabs(*(matA+i+ *(P+i) * lda)))
                *(P+i)=j;
        }
        if (fabs(*(matA+i+ *(P+i) * lda))<1e-5) {
            printf("warning: matrix is singular in row %d (function: dense_LU_small) \n", *(P+i));
            *(P+i)= -1;
            continue;
        }
        else if (*(P+i)!=i) {
            for (j=0;j<n;++j) {
                temp = *(matA+i*lda+j);
                *(matA+i*lda+j)=*(matA + *(P+i) * lda + j);
                *(matA + *(P+i) * lda + j)=temp;
            }
        }
        *length = m-i-1;
        drscl( length, matA+i*lda+i,matA+i+(i+1)*lda, &lda);
        cblas_dger(CblasRowMajor,m-i-1,n-i-1,-1.0,matA+(i+1)*lda+i,lda,matA+i*lda+i+1,1,matA+(i+1)*lda+i+1,lda);
    }
    free(length);
}

void dense_LU_block(int n,int b, dense_matrix matA, int *P)
{
    int i,j,k,end, *length;
    double temp;
    length=malloc(sizeof(int));
    for (i=0;i<(n-1);i+=b) {
        end=i+b;
        if (end>n)
            end=n;
        dense_LU_small(end-i,n-i,n,matA+i*n+i, P+i);
        for (k=i;k<end;++k) {
            *(P+k) += i;
            if (*(P+k) != k) {
                for (j=0;j<i;++j) {
                    temp = *(matA+k*n+j);
                    *(matA+k*n+j)=*(matA + *(P+k) * n + j);
                    *(matA + *(P+k) * n + j)=temp;
                }
                for (j=end;j<n;++j) {
                    temp = *(matA+k*n+j);
                    *(matA+k*n+j)=*(matA + *(P+k) * n + j);
                    *(matA + *(P+k) * n + j)=temp;
                }
            }
        }
        *length=n-end;
        cblas_dtrsm(CblasRowMajor,CblasLeft,CblasLower,CblasNoTrans,CblasUnit,end-i,n-end,1,matA+i*n+i,n,matA+i*n+end,n);
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,n-end,n-end,end-i,-1,matA+end*n+i,n,matA+i*n+end,n,1,matA+end*n+end,n);
    }
    free(length);
}

void sparse_row_switch(int n, struct sparse_CRS *matA, int row1,int row2)
{
    int i,j,begin1,begin2,length1,length2,*tempcol;
    double *tempval;
    if (row2<row1) {
        printf("row2 must be larger than row1 (function: sparse_row_switch)");
        return;
    }
    if (row1==row2)
        return;
    begin1 = *(matA->rowptr+row1);
    begin2 = *(matA->rowptr+row2);
    length1 = *(matA->rowptr+row1+1)-*(matA->rowptr+row1);
    length2 = *(matA->rowptr+row2+1)-*(matA->rowptr+row2);
    tempcol=malloc(sizeof(int));
    tempval=malloc(sizeof(double));
    for (i=0;i<length1 && i<length2;++i) {
        *tempcol= *(matA->colnumber+begin1+i);
        *(matA->colnumber+begin1+i)=*(matA->colnumber+begin2+i);
        *(matA->colnumber+begin2+i)=*tempcol;
        *tempval= *(matA->value+begin1+i);
        *(matA->value+begin1+i)=*(matA->value+begin2+i);
        *(matA->value+begin2+i)=*tempval;
    }
    if (length1>length2) {
        tempcol=realloc(tempcol,(length1-length2)*sizeof(int));
        tempval=realloc(tempval,(length1-length2)*sizeof(double));
        for (i=length2;i<length1;++i) {
            *(tempcol+i-length2)=*(matA->colnumber+begin1+i);
            *(tempval+i-length2)=*(matA->value+begin1+i);
        }
        for (i=begin1+length1;i<begin2+length2;++i) {
            *(matA->colnumber+i-length1+length2)=*(matA->colnumber+i);
            *(matA->value+i-length1+length2)=*(matA->value+i);
        }
        for (i=length2;i<length1;++i) {
            *(matA->colnumber+begin2-length1+length2+i)=*(tempcol+i-length2);
            *(matA->value+begin2-length1+length2+i)=*(tempval+i-length2);
        }
    }
    else if (length1<length2) {
        tempcol=realloc(tempcol,(length2-length1)*sizeof(int));
        tempval=realloc(tempval,(length2-length1)*sizeof(double));
        for (i=length1;i<length2;++i) {
            *(tempcol+i-length1)=*(matA->colnumber+begin2+i);
            *(tempval+i-length1)=*(matA->value+begin2+i);
        }
        for (i=begin2+length1-1;i>=begin1+length1;--i) {
            *(matA->colnumber+i-length1+length2)=*(matA->colnumber+i);
            *(matA->value+i-length1+length2)=*(matA->value+i);
        }
        for (i=length1;i<length2;++i) {
            *(matA->colnumber+begin1+i)=*(tempcol+i-length1);
            *(matA->value+begin1+i)=*(tempval+i-length1);
        }
    }
    for (i=row1+1;i<=row2;++i) {
        *(matA->rowptr+i) = *(matA->rowptr+i) - length1 + length2;
    }
    free(tempcol);
    free(tempval);
}

void sparse_LU(int n,struct sparse_CRS *matA, int *P)
{
    int i,j,k,firstj,firstcount, row,rowcol,count,temp,colj;
    char *number;
    number=malloc(8*sizeof(char));
    for (i=0;i<n;++i) {
        *(P+i)=i;
    }
    for (i=0;i<n;++i) {
        count=i;
        for (j=i;j<n;++j) {
            firstj=*(matA->rowptr+j);
            while (firstj<*(matA->rowptr+j+1) && *(matA->colnumber+firstj)<i)
                ++firstj;
            if (firstj<*(matA->rowptr+j+1) && *(matA->colnumber+firstj)==i) {
                sparse_row_switch(n,matA,count,j);
                temp=*(P+j);
                *(P+j)=*(P+count);
                *(P+count)=temp;
                firstcount=*(matA->rowptr+count);
                while (*(matA->colnumber+firstcount)<i)
                    ++firstcount;
                row=*(matA->rowptr+i);
                while (*(matA->colnumber+row)<i)
                    ++row;
                if (fabs(*(matA->value+row))<fabs(*(matA->value+firstcount))) {
                    sparse_row_switch(n,matA,i,count);
                    temp=*(P+i);
                    *(P+i)=*(P+count);
                    *(P+count)=temp;
                    row=*(matA->rowptr+i);
                    while (*matA->colnumber+row<i)
                        ++row;
                }
                ++count;
            }
        }
        if (count==i) {
            printf("Sparse matrix is singular (function: sparse_LU)\n");
            return;
        }
        for (j=i+1;j<count;++j) {
            firstj=*(matA->rowptr+j);
            while (*(matA->colnumber+firstj) < i)
                ++firstj;
            *(matA->value+firstj) /= *(matA->value+row);
            rowcol=row;
            while (rowcol+1 < *(matA->rowptr+i+1)) {
                ++rowcol;
                colj=firstj;
                while (*(matA->colnumber+colj)<*(matA->colnumber+rowcol) && colj<*(matA->rowptr+j+1))
                    ++colj;
                if (colj<*(matA->rowptr+j+1) && *(matA->colnumber+colj)==*(matA->colnumber+rowcol))
                    *(matA->value+colj) = *(matA->value+colj) - *(matA->value+firstj) * *(matA->value+rowcol);
                else {
                    for (k=j+1;k<(n+1);++k)
                        *(matA->rowptr+k) +=1;
                    matA->colnumber=realloc(matA->colnumber,*(matA->rowptr+n) * sizeof(int));
                    matA->value=realloc(matA->value,*(matA->rowptr+n) * sizeof(double));
                    for (k=*(matA->rowptr+n)-1;k>colj;--k) {
                        *(matA->colnumber+k)=*(matA->colnumber+k-1);
                        *(matA->value+k)=*(matA->value+k-1);
                    }
                    *(matA->colnumber+colj)=*(matA->colnumber+rowcol);
                    *(matA->value+colj)=*(matA->value+rowcol) * *(matA->value+firstj) * (-1.0);
                }
            }
        }
    }
    /*  for(i=0;i<n;++i){
        if(*(P+i)<i)
          *(P+i)=*(P+*(P+i));
      }*/
}
