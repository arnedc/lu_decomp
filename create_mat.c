#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"

extern long seed;

/**
 * @brief creates a random dense matrix of size n x n. The values of the dense matrix are in the area [-max, +max[.
 *
 * @param seed this is the seed used for the random distribution of the values of the matrix
 * @param n defines the size of the dense matrix (n x n)
 * @param max defines the interval in which the values of the matrix can lie
 * @return dense_matrix
 **/
dense_matrix create_matrix(long seed, int n, double max)
{
    int i,j;
    srand48(seed);
    double * matrix;
    matrix = (double *) malloc(n * n * sizeof(double));
    if (matrix==NULL) {
        printf("Unable to create dense matrix (function: create_matrix)");
        return NULL;
    }
    for (i=0; i<n;++i) {
        for (j=0;j<n;++j) {
            *(matrix + j +i*n) = drand48()*2.0*max-max;
        }
    }
    return matrix;
}

void destroy_dense(void *matrix)
{
    free(matrix);
    return;
}

dense_matrix create_symm(long seed, int n, double max)
{
    int i,j;
  srand48(seed);
  double * matrix;
  matrix = (double *) calloc((n*n+n)/2, sizeof(double));
  if (matrix==NULL){
    printf("Unable to create dense matrix (function: create_symm)");
    return NULL;
  }
  for (i=0; i<n;++i) {
    for(j=0;j<=i;++j) {
      *(matrix + j + (i*i+i)/2) = drand48()*2.0*max-max;
    }
  }
  return matrix;
}

dense_matrix create_vect(long seed, int n, double max)
{
    int i;
    srand48(seed);
    double * matrix;
    matrix = (double *) malloc(n * sizeof(double));
    if (matrix==NULL) {
        printf("Unable to create dense matrix (function: create_matrix)");
        return NULL;
    }
    for (i=0; i<n;++i) {
       
            *(matrix + i) = drand48()*2.0*max-max;
       
    }
    return matrix;
}

dense_matrix CRS2dense(int n, struct sparse_CRS *matA)
{
    int i,j,col;
    double *matrix;
    matrix = (double *) malloc(n * n * sizeof(double));
    if (matrix==NULL) {
        printf("Unable to create dense matrix (function: create_matrix)");
        return NULL;
    }
    for (i=0;i<n;++i) {
        col=*(matA->rowptr+i);
        for (j=0;j<n;++j) {
            if (col < *(matA->rowptr+i+1) && *(matA->colnumber+col)==j) {
                *(matrix+i*n+j)=*(matA->value+col);
                ++col;
            }
            else
                *(matrix+i*n+j)=0;
        }
    }
    return matrix;
}

struct sparse_CRS *create_sparse_CRS(long seed, int n, double fill, double max)
{
    int i,j, *c, nz,pos,rp;
    double *v;
    srand48(seed);
    struct sparse_CRS *mat;
    mat=(struct sparse_CRS *) malloc(sizeof (struct sparse_CRS));
    if (mat==NULL) {
        printf("Unable to create sparse matrix (function: create_sparse_CRS)");
        return NULL;
    }
    mat->rowptr=calloc(n+1,sizeof(int));
    if (mat->rowptr==NULL) {
        printf("Unable to create sparse matrix (function: create_sparse_CRS)");
        return NULL;
    }
    mat->colnumber=malloc(0);
    mat->value=malloc(0);
    *(mat->rowptr) = 0;
    rp=0;
    for (i=0; i<n; ++i) {
        nz=(int) (drand48()*n*fill/50);
        if (nz==0)
            nz=1;
        rp+=nz;
        *(mat->rowptr + i + 1) = rp;
        mat->colnumber=(int *) realloc(mat->colnumber,rp*sizeof(int));
        if (mat->colnumber==NULL) {
            printf("Unable to create sparse matrix (function: create_sparse_CRS)");
            return NULL;
        }
        mat->value=(double *) realloc(mat->value,rp*sizeof(double));
        if (mat->value==NULL) {
            printf("Unable to create sparse matrix (function: create_sparse_CRS)");
            return NULL;
        }
        pos=0;
        c=mat->colnumber + *(mat->rowptr+i);
        v=mat->value + *(mat->rowptr+i);
        for (j=0; j<nz; ++j) {
            pos += j==0 ? (int) (drand48()*(n-nz)) : (int) (1+drand48()*(n-nz+j-pos));
            *c++ = pos;
            *v++ = drand48()*2.0*max-max;
        }
    }
    return mat;
}

void destroy_sparse_CRS(struct sparse_CRS *matrix)
{
    free(matrix->colnumber);
    free(matrix->rowptr);
    free(matrix->value);
    free(matrix);
    return;
}

/**
 * @brief creates a random sparse matrix of size n x n. The values of the sparse matrix lie in the interval [-max, +max[.
 *
 * The sparse matrix is a pointer to a row structure, incrementing the sparse_matrix makes it point to the next row.
 * A row is defined as a structure that consists of a pointer to colnumber (integers), which holds the indices of the non-zero elements
 * and a pointer to values (doubles), which holds the values of the matrix at the corresponding column index.
 *
 * @param seed this is the seed used for the random distribution of the values of the matrix
 * @param n defines the size of the dense matrix (n x n)
 * @param max defines the interval in which the values of the matrix can lie
 * @return sparse_matrix
 **/
sparse_matrix create_sparse_matrix(long seed, int n, double max)
{
    int i,j, *c, nz,pos;
    double *v;
    srand48(seed);
    struct row *p;
    p = (struct row *) calloc(n, sizeof(struct row));
    if (p==NULL) {
        printf("Unable to create sparse matrix (function: create_sparse_matrix)");
        return NULL;
    }
    for (i=0; i<n;++i) {
        nz=(int) (drand48()*n);
        c=(int *) calloc(nz+1,sizeof(int));
        if (c==NULL) {
            printf("Unable to create sparse matrix (function: create_sparse_matrix)");
            return NULL;
        }
        v=(double *) calloc(nz+1,sizeof(double));
        if (v==NULL) {
            printf("Unable to create sparse matrix (function: create_sparse_matrix)");
            return NULL;
        }
        (p+i)->colnumber = c;
        (p+i)->value = v;
        pos=0;
        for (j=0;j<nz;++j) {
            pos += j==0 ? (int) (drand48()*(n-nz)) : (int) (1+drand48()*(n-nz+j-pos));
            *c++ = pos;
            *v++ = drand48()*2.0*max-max;
        }
        *c=-1;
        *v=-1;
    }
    return p;
}

void destroy_sparse(int n, sparse_matrix matrix)
{
    int i;
    for (i=0;i<n;++i) {
        free((matrix+i)->colnumber);
        free((matrix+i)->value);
    }
    free(matrix);
    return;
}


/**
 * @brief creates a random sparse matrix of size n x n. The values of the sparse matrix lie in the interval [-max, +max[.
 *
 * The sparse matrix is a pointer to a cellrow structure, incrementing the sparse_matrix makes it point to the next cellrow.
 * A cellrow is defined as a structure that consists of a pointer to a cell. Every row starts with a head, which is an empty cell.
 * A cell is defined as a structure holding the column index (colnumber), the value of the cell (value) and a pointer to the next cell.
 * If the pointer to the next cell is a NULL-pointer, there are no more non-zero elements in the row.
 *
 * @param seed this is the seed used for the random distribution of the values of the matrix
 * @param n defines the size of the dense matrix (n x n)
 * @param max defines the interval in which the values of the matrix can lie
 * @return sparse2_matrix
 **/
sparse2_matrix create_sparse2_matrix(long seed, int n, double max)
{
    int i,j,nz,pos, count=0;
    struct cell *huidig;
    srand48(seed);
    struct cellrow *mat;
    mat = (struct cellrow *) calloc(n+1,sizeof(struct cellrow));
    if (mat==NULL) {
        printf("Unable to create sparse matrix (function: create_sparse2_matrix)");
        return NULL;
    }
    for (i=0; i<n;++i) {
        nz=(int) (drand48()*n);
        pos=0;
        count +=nz;
        (mat+i)->head=malloc(sizeof (struct cell));
        if ((mat+i)->head==NULL) {
            printf("Unable to create sparse matrix (function: create_sparse2_matrix)");
            return NULL;
        }
        huidig=(mat+i)->head;
        huidig->colnumber=-1;
        huidig->value=0;
        for (j=0;j<nz;++j) {
            pos += j==0 ? (int) (drand48()*(n-nz)) : (int) (1+drand48()*(n-nz+j-pos));
            huidig->next=malloc(sizeof (struct cell));
            if (huidig->next==NULL) {
                printf("Unable to create sparse matrix (function: create_sparse2_matrix)");
                return NULL;
            }
            huidig=huidig->next;
            huidig->colnumber = pos;
            huidig->value = drand48()*2.0*max-max;
        }
        huidig->next=NULL;
    }
    printf("Sparse matrix created %g %% filled\n",(double) (count)/n/n*100);
    return mat;
}

void destroy_sparse2(int n, sparse2_matrix matrix)
{
    int i;
    struct cell *huidig;
    for (i=0;i<n;++i) {
        while ((matrix+i)->head != NULL) {
            huidig=(matrix+i)->head->next;
            free((matrix+i)->head);
            (matrix+i)->head=huidig;
        }
        free((matrix+i)->head);
    }
    free(matrix);
    return;
}

dense_matrix sparse2dense(int n, sparse2_matrix sparse)
{
    struct cell *huidig;
    double *matrix;
    matrix = (double *) malloc(n * n * sizeof(double));
    if (matrix==NULL) {
        printf("Unable to create dense matrix (function: sparse2dense)");
        return NULL;
    }
    int i,j;
    for (i=0;i<n;++i) {
        huidig=(sparse+i)->head->next;
        for (j=0;j<n;++j) {
            if (huidig != NULL && j== huidig->colnumber) {
                *(matrix + j +i*n)= huidig->value;
                huidig=huidig->next;
            }
            else
                *(matrix + j +i*n)=0;
        }
    }
    return matrix;
}

sparse_matrix sparse2sparse(int n, sparse2_matrix sparse)
{
    struct cell *huidig;
    struct row *matrix;
    matrix = (struct row *) calloc(n+1, sizeof(struct row));
    if (matrix==NULL) {
        printf("Unable to create sparse matrix (function: sparse2sparse)");
        return NULL;
    }
    int i,j,k,*c;
    double *v;
    for (i=0;i<n;++i) {
        huidig=(sparse+i)->head->next;
        for (j=0;huidig!=NULL;++j)
            huidig=huidig->next;
        huidig=(sparse+i)->head->next;
        c=(int *) calloc(j+1,sizeof(int));
        if (c==NULL) {
            printf("Unable to create sparse matrix (function: sparse2sparse)");
            return NULL;
        }
        v=(double *) calloc(j+1,sizeof(double));
        if (v==NULL) {
            printf("Unable to create sparse matrix (function: sparse2sparse)");
            return NULL;
        }
        (matrix+i)->colnumber = c;
        (matrix+i)->value = v;
        for (k=0;k<j;++k) {
            *c = huidig->colnumber;
            c++;
            *v = huidig->value;
            v++;
            huidig=huidig->next;
        }
        *c=-1;
        *v=-1;
    }
    return matrix;
}

dense_matrix symm2dense(int n, double *symm)
{ 
  int i,j;
  double * matrix;
  matrix = (double *) calloc(n*n, sizeof(double));
  if (matrix==NULL){
    printf("Unable to create dense matrix (function: symm2dense)");
    return NULL;
  }
  for (i=0; i<n;++i) {
    for(j=0;j<=i;++j) {
      *(matrix+i*n+j)=*(matrix+j*n+i)= *(symm + j + (i*i+i)/2);
    }
  }
  return matrix;
}

dense_matrix delete_row_symm(int n, int row, dense_matrix matA)
{
  int i,j;
  for (i=row+1;i<n;++i){
    for (j=0;j<row;++j){
      *(matA+(i*i-i)/2+j)=*(matA+(i*i+i)/2+j);
    }
    for (j=row+1;j<=i;++j){
      *(matA+(i*i-i)/2+j-1)=*(matA+(i*i+i)/2+j);
    }
  }
  matA=realloc(matA,(n*n-n)/2*sizeof(double));
  return matA;
}

dense_matrix symm2RFP(int n, double *symm)
{
  int i, j;
  double *RFP;
  RFP=(double *) calloc((n*n+n)/2, sizeof(double));
  if(RFP==NULL){
    printf("unable to create RFP matrix (function: symm2RFP) \n");
    return;
  }
  if(n % 2 ==0){
    for (i=0;i<=n/2;++i){
      for (j=0;j<n/2;++j){
	if(j<i)
	  *(RFP+i*n/2+j)=*(symm+j+(i*i-i)/2);
	else
	  *(RFP+i*n/2+j)=*(symm+((j+n/2)*(j+n/2+1))/2+i+n/2);
      }
    }
    for(i=n/2+1;i<n+1;++i){
      for (j=0;j<n/2;++j){
	*(RFP+i*n/2+j)=*(symm+(i*i-i)/2+j);
      }
    }
  }
  else{
    for (i=0;i<(n+1)/2;++i){
      for (j=0;j<(n+1)/2;++j){
	if(j>i)
	  *(RFP+i*(n+1)/2+j)=*(symm+((j+(n+1)/2-1)*(j+(n+1)/2))/2+i+(n+1)/2);
	else
	  *(RFP+i*(n+1)/2+j)=*(symm+j+(i*i+i)/2);
      }
    }
    for(i=(n+1)/2;i<n;++i){
      for (j=0;j<(n+1)/2;++j){
	*(RFP+i*(n+1)/2+j)=*(symm+(i*i+i)/2+j);
      }
    }
  }
  return RFP;
}
