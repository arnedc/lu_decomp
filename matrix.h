long seed;

struct cell {
  int colnumber;
  double value;
  struct cell *next;
};

struct cellrow {
  struct cell *head;
};

struct row {
  int *colnumber;
  double *value;
};

struct sparse_CRS {
  int *colnumber;
  double *value;
  int *rowptr;
};

typedef struct cellrow *sparse2_matrix;
typedef double *dense_matrix, *dense_tri;
typedef struct row *sparse_matrix;

dense_matrix create_matrix(long,int, double);
void destroy_dense(void *);
dense_matrix multiply_matrix(int, double *, double *);
dense_matrix multiply_blas_matrix(int, double *, double *);
dense_matrix create_vect(long seed, int n, double max);

sparse_matrix create_sparse_matrix(long, int, double);
void destroy_sparse(int,sparse_matrix);
sparse_matrix multiply_sparse_matrix(int, struct row *, struct row *);
struct sparse_CRS *create_sparse_CRS(long, int, double, double);
void destroy_sparse_CRS(struct sparse_CRS *);
struct sparse_CRS *multiply_CRS_dense(int, struct sparse_CRS *, struct sparse_CRS *);
struct sparse_CRS *multiply_CRS_sparse(int, struct sparse_CRS *, struct sparse_CRS *);

sparse2_matrix create_sparse2_matrix(long, int, double);
void destroy_sparse2(int,sparse2_matrix);
sparse2_matrix multiply_sparse2_sparse(int , struct cellrow *, struct cellrow *);
sparse2_matrix multiply_sparse2_dense(int, sparse2_matrix, sparse2_matrix);

dense_matrix sparse2dense(int,sparse2_matrix);
sparse_matrix sparse2sparse(int,sparse2_matrix);
struct sparse_CRS *sparse2CRS(int, sparse2_matrix);
dense_matrix symm2RFP(int n, double *symm);

void printdense(int, dense_matrix, char*);
void printsparse(int, sparse_matrix, char *);
void printsparse2(int, sparse2_matrix, char *);
void printsparse_CRS(int, struct sparse_CRS *, char *);
void printvect(int n, double *perm, char *filename);
void printsymm(int, dense_matrix, char *);
void printperm(int, int *, char *);
void printsymm_RFP(int n, dense_matrix mat, char *filename);

dense_matrix create_symm(long, int, double);
dense_matrix symm2dense(int, double *);

dense_matrix multiply_symm_matrix(int, dense_matrix, dense_matrix);

void transpose(int,dense_matrix);

void dense_LU(int,dense_matrix, dense_tri, dense_tri, int *);
int lapack_LU(int,dense_matrix,int *);
int lapack_symm_LU(int n,dense_matrix matA, int *P);
void dense_LU_small(int,int,int,dense_matrix, int *);
void dense_LU_block(int,int,dense_matrix, int *);
void dense_cholesky(int n,dense_matrix matA);
dense_matrix delete_row_symm(int n, int row, dense_matrix matA);
void dense_RFP_cholesky(int n,dense_matrix matA);
void dense_RFP_solve(int n, dense_matrix matA, double *matB);

void sparse_row_switch(int n, struct sparse_CRS *matA, int row1,int row2);
void sparse_LU(int n,struct sparse_CRS *matA, int *P);

dense_matrix CRS2dense(int n, struct sparse_CRS *matA);

int import_symm(char *filename,dense_matrix *symm);
int import_vect(char *filename,double **vect);

int lapack_symm_solve(int n,dense_matrix matA,int *P, double *matB);
double *dense_solve(int n, dense_matrix matA, int *P, double *matB);
double *dense_cholsolve(int n, dense_matrix matA, double *matB);

