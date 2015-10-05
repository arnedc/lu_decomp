#include <stdio.h>
#include "matrix.h"
#include <stdlib.h>

int import_symm(char *filename,dense_matrix *symm)
{
    FILE *fd;
    char *line;
    int row, col,old;
    int value;
    fd = fopen(filename,"r");
    if (fd==NULL)
        printf("error reading file");
    row=col=0;
    value=0.0;
    line=malloc(80*sizeof(char));
    while (fgets(line, 80, fd) != NULL) {
        if (sscanf(line,"%d %d %d",&row,&col,&value)!=3)
            printf("wrong input");
    }
    fclose(fd);
    *symm=calloc((row*row+row)/2,sizeof(double));
    fd = fopen(filename,"r");
    if (fd==NULL)
        printf("error reading file");
    old=0;
    while (fgets(line, 80, fd) != NULL) {
        if (sscanf(line,"%d %d %d",&row,&col,&value)!=3)
            printf("wrong input");
        old=(row*row-row)/2+col;
        *(*symm+old-1)=(double) value;
    }
    fclose(fd);
    free(line);
    return row;
}

int import_vect(char *filename,double **vect)
{
    FILE *fd;
    char *line;
    line=malloc(80*sizeof(char));
    int count=0;
    double value=0.0;
    fd = fopen(filename,"r");
    if (fd==NULL)
        printf("error reading file");
    while (fgets(line, 80, fd) != NULL) {
        ++count;
        if (sscanf(line,"%lf",&value)!=1)
            printf("wrong input");
    }
    fclose(fd);
    *vect=calloc(count,sizeof(double));
    fd = fopen(filename,"r");
    if (fd==NULL)
        printf("error reading file");
    count=0;
    while (fgets(line, 80, fd) != NULL) {
        if (sscanf(line,"%lf",&value)!=1)
            printf("wrong input");
        *(*vect+count)=value;
        ++count;
    }
    fclose(fd);
    free(line);
    return count;
}
