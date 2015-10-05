#include <stdio.h>
#include "matrix.h"

void printdense(int n, dense_matrix mat, char *filename)
{
    FILE *fd;
    fd = fopen(filename,"w");
    if (fd==NULL)
        printf("error creating file");
    int i,j;
    for (i=0;i<n;++i) {
        fprintf(fd,"[\t");
        for (j=0;j<n;++j) {
            fprintf(fd,"%7.3f\t",*mat++);
        }
        fprintf(fd,"]\n");
    }
    fclose(fd);
}

void printsymm(int n, dense_matrix mat, char *filename)
{
    FILE *fd;
    fd = fopen(filename,"w");
    if (fd==NULL)
        printf("error creating file");
    int i,j;
    for (i=0;i<n;++i) {
        fprintf(fd,"[\t");
        for (j=0;j<n;++j) {
            if (j>i)
                fprintf(fd,"%7.3f\t",*(mat+i+(j*j+j)/2));
            else
                fprintf(fd,"%7.3f\t",*(mat+j+(i*i+i)/2));
        }
        fprintf(fd,"]\n");
    }
    fclose(fd);
}

void printperm(int n, int *perm, char *filename)
{
    FILE *fd;
    fd = fopen(filename,"w");
    if (fd==NULL) {
        printf("error creating file");
        return;
    }
    int i,j;
    fprintf(fd,"[\t");
    for (i=0;i<n;++i) {
        fprintf(fd,"%d\t",*perm++);
    }
    fprintf(fd,"]\n");
    fclose(fd);
}

void printvect(int n, double *perm, char *filename)
{
    FILE *fd;
    fd = fopen(filename,"w");
    if (fd==NULL) {
        printf("error creating file");
        return;
    }
    int i;
    for (i=0;i<n;++i) {
        fprintf(fd,"%10.6f\n",*perm++);
    }
    fclose(fd);
}

void printsparse_CRS(int n, struct sparse_CRS *mat, char *filename)
{
    FILE *fd;
    int max,count;
    fd = fopen(filename,"w");
    if (fd==NULL)
        printf("error creating file");
    int i,j;
    count=0;
    for (i=0;i<n;++i) {
        fprintf(fd,"[\t");
        count= *(mat->rowptr+i);
        max= *(mat->rowptr+i+1);
        for (j=0;j<n;++j) {
            if (count<max && j== *(mat->colnumber+count)) {
                fprintf(fd,"%7.3f\t",*(mat->value+count));
                ++count;
            }
            else
                fprintf(fd,"%7.3f\t",0.0);
        }
        fprintf(fd,"]\n");
    }
    fclose(fd);
}

void printsparse(int n, sparse_matrix mat, char *filename)
{
    FILE *fd;
    struct row huidig;
    fd = fopen(filename,"w");
    if (fd==NULL)
        printf("error creating file");
    int i,j;
    for (i=0;i<n;++i) {
        fprintf(fd,"[\t");
        huidig=*(mat+i);
        for (j=0;j<n;++j) {
            if (j== *huidig.colnumber) {
                fprintf(fd,"%7.3f\t",*huidig.value);
                ++(huidig.value);
                ++(huidig.colnumber);
            }
            else
                fprintf(fd,"%7.3f\t",0.0);
        }
        fprintf(fd,"]\n");
    }
    fclose(fd);
}

void printsparse2(int n, sparse2_matrix mat, char *filename)
{
    FILE *fd;
    struct cell *huidig;
    fd = fopen(filename,"w");
    if (fd==NULL)
        printf("error creating file");
    int i,j;
    for (i=0;i<n;++i) {
        fprintf(fd,"[\t");
        huidig=(mat+i)->head->next;
        for (j=0;j<n;++j) {
            if (huidig != NULL && j== huidig->colnumber) {
                fprintf(fd,"%7.3f\t",huidig->value);
                huidig=huidig->next;
            }
            else
                fprintf(fd,"%7.3f\t",0.0);
        }
        fprintf(fd,"]\n");
    }
    fclose(fd);
}

void printsymm_RFP(int n, dense_matrix mat, char *filename)
{
  FILE *fd;
  fd = fopen(filename,"w");
  if (fd==NULL)
    printf("error creating file");
  int i,j,rows,cols;
  if(n % 2 == 0){
    rows=n+1;
    cols=n/2;
  }
  if(n % 2 == 1){
    rows=n;
    cols=(n+1)/2;
  }
  for (i=0;i<rows;++i) {
    fprintf(fd,"[\t");
    for (j=0;j<cols;++j) {
      fprintf(fd,"%7.3f\t",*mat++);
    }
    fprintf(fd,"]\n");
  }
  fclose(fd);
}
