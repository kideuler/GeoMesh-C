#include "GeoMesh.h"

// creating matrices
struct DoubleMatrix DoubleMatrix_create(int nrows, int ncols){
    struct DoubleMatrix M;
    M.nrows = nrows;
    M.ncols = ncols;
    M.data = (double*) malloc(nrows * ncols * sizeof(double));
    return M;
}
struct IntMatrix IntMatrix_create(int nrows, int ncols){
    struct IntMatrix M;
    M.nrows = nrows;
    M.ncols = ncols;
    M.data = (int*) malloc(nrows * ncols * sizeof(int));
    return M;
}

// printing matrices
void DoubleMatrix_print(struct DoubleMatrix* M){
    for (int i = 0; i<M->nrows; i++){
        for (int j = 0; j<M->ncols; j++){
            printf("%f ",M->data[i*M->ncols + j]);
        }
        printf("\n");
    }
}
void IntMatrix_print(struct IntMatrix* M){
    for (int i = 0; i<M->nrows; i++){
        for (int j = 0; j<M->ncols; j++){
            printf("%d ",M->data[i*M->ncols + j]);
        }
        printf("\n");
    }
}


double DoubleMatrix_min(struct DoubleMatrix* M, int col){
    double val = M->data[col];
    for (int i = 1; i<M->nrows; i++){
        if (M->data[i*M->ncols + col] < val){
            val = M->data[i*M->ncols + col];
        }
    }
    return val;
}
double DoubleMatrix_max(struct DoubleMatrix* M, int col){
    double val = M->data[col];
    for (int i = 1; i<M->nrows; i++){
        if (M->data[i*M->ncols + col] > val){
            val = M->data[i*M->ncols + col];
        }
    }
    return val;
}

double drand(double low, double high){
  double d;
  d = (double) rand() / ((double) RAND_MAX + 1);
  return (low + d * (high - low));
}
struct DoubleMatrix DoubleMatrix_create_Random(int nrows, int ncols, double low, double high){
    struct DoubleMatrix M = DoubleMatrix_create(nrows,ncols);
    for (int i = 0; i<nrows; i++){
        for (int j = 0; j<ncols; j++){
            M.data[i*M.ncols + j] = drand(low,high);
        }
    }
    return M;
}