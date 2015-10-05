#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"
#include <time.h>
#include <mkl_cblas.h>
#include <mkl_lapacke.h>
#include <string.h>

int main(int argc, char **argv) {
    //struct sparse_CRS *matA;
    double *RFPA, *symmA, **ptrA, *y, **ptry, *sol, *fullA;
    double *matA, *matcopy;
    int i,m,n, *P, *Pdense, info,count, b;
    clock_t c0, c1,c2;
    double max=1.0, fill=50;
    seed=time(NULL);
    ptrA=&symmA;
    ptry=&y;
    char preamble[100];
    if (argc < 4){
      printf("Invalid number of input parameters. \n");
      printf("./LU_decomp <matrix size> <block size> <preamble>\n");
      return;
    }
    n=atoi(*++argv);
    printf("Matrix size: %d\n",n);
    b=atoi(*++argv);
    printf("Block size: %d \n",b);
    strcpy(preamble, *++argv);
    printf("Preamble: %s \n",preamble);
    
   /*    n=import_symm("/home/arnedc/projects/3SNP_WRW.txt",ptrA);
        m=import_vect("/home/arnedc/projects/3SNP_WRy.txt",ptry);
        symmA=create_symm(seed,n,2.0);
        y=create_vect(seed,n,2.0);
        printsymm(n,symmA,"Asymm");
        printvect(n,y,"ymat");
        RFPA=symm2RFP(n,symmA);
        printsymm_RFP(n,RFPA,"Asymmcop");
        dense_cholesky(n,symmA);
        fullA=dense_cholsolve(n,symmA,y);
        dense_RFP_cholesky(n,RFPA);
        dense_RFP_solve(n,RFPA,y);
        printsymm_RFP(n,RFPA,"ARFPchol");
        printsymm(n,symmA,"Asymmchol");
        printvect(n,fullA,"chol_sol");
        printvect(n,y,"RFP_sol");*/

    /* matB=calloc(n*n,sizeof(double));
      if(matB==NULL){
        printf("Unable to create dense triangular matrix (function: dense_LU)");
        return;
      }
      cblas_dcopy(n*n,matA,1,matB,1);
      upA=calloc(n*n,sizeof(double));
      if(upA==NULL){
        printf("Unable to create dense triangular matrix (function: dense_LU)");
        return;
      }*/

    matcopy= (double *) calloc(n*n, sizeof(double));
    matA=create_matrix(seed,n,5.0);
    cblas_dcopy(n*n,matA,1,matcopy,1);

    P=(int *) calloc(n,sizeof(int));
    if (P==NULL) {
        printf("Unable to create dense triangular matrix (function: dense_LU)");
        return;
    }
    Pdense=(int *) calloc(n,sizeof(int));
    if (Pdense==NULL) {
        printf("Unable to create dense triangular matrix (function: dense_LU)");
        return;
    }
    c0=clock();
    dense_LU_small(n,n,n,matA,Pdense);
    //lapack_LU(n,matB,Pdense);
    c1=clock();
    printf("CPU time for dense LU:           %.8f \n",(double) (c1-c0)/CLOCKS_PER_SEC);
    char filename[100];
    filename[0]='\0';
    strcat(filename,preamble);
    strcat(filename,"_LU.txt");
    //printf("Filename unblocked: %s\n", filename);
    printdense(n,matA,filename);
    printperm(n,Pdense,"permLUA.txt");
    cblas_dcopy(n*n,matcopy,1,matA,1);
    c0=clock();
    dense_LU_block(n,b,matA,P);
    c2=clock();
    printf("CPU time for blockwise dense LU:           %.8f \n",(double) (c2-c0)/CLOCKS_PER_SEC);
    lapack_LU(n,matcopy,Pdense);
    c1=clock();

    printf("CPU time for dense LU (lapack):   %.8f \n",(double) (c1-c2)/CLOCKS_PER_SEC);

    filename[0]='\0';
    strcat(filename,preamble);
    strcat(filename,"_BlockLU.txt");
    printdense(n,matA,filename);
    filename[0]='\0';
    strcat(filename,preamble);
    strcat(filename,"_LAPACKLU.txt");
    printdense(n,matcopy,filename);
    printperm(n,P,"permblock.txt");
    printperm(n,Pdense,"LAPACKperm.txt");
    destroy_dense(P);
    destroy_dense(Pdense);

    destroy_dense(matA);
    destroy_dense(matcopy);

    return 0;
}
