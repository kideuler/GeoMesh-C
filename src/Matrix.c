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

// resizing functions
void DoubleMatrix_resize(struct DoubleMatrix* M, int newsz){
    double* data = (double*)malloc(newsz * M->ncols * sizeof(double));
    int size = (newsz>M->nrows)?M->nrows:newsz;
    for (int i = 0; i< size; i++){
        for (int j = 0; j<M->ncols; j++){
            data[M->ncols*i + j] = M->data[M->ncols*i + j];
        }
    }
    M->nrows = newsz;
    free(M->data);
    M->data = data;
}
void IntMatrix_resize(struct IntMatrix* M, int newsz){
    int* data = (int*)malloc(newsz * M->ncols * sizeof(int));
    int size = (newsz>M->nrows)?M->nrows:newsz;
    for (int i = 0; i<size; i++){
        for (int j = 0; j<M->ncols; j++){
            data[M->ncols*i + j] = M->data[M->ncols*i + j];
        }
    }
    M->nrows = newsz;
    free(M->data);
    M->data = data;
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

void load_lakeSuperior(struct IntMatrix* segments, struct DoubleMatrix* coords){
    FILE *fid;
    int nv,nsegs,v1,v2,v3;
    float x, y, z;
    free(segments->data);
    free(coords->data);

    fid = fopen("data/lake.dat","r");
    fscanf(fid,"%d",&nv);
    *coords = DoubleMatrix_create(nv,2);
    for(int i=0; i<nv; i++){
        fscanf(fid,"%f %f %f", &x,&y,&z);
        coords->data[2*i] = x;
        coords->data[2*i+1] = y;
    }

    fscanf(fid,"%d",&nsegs);
    *segments = IntMatrix_create(nsegs,2);
    for(int i=0; i<nsegs; i++){
        fscanf(fid,"%d %d %d", &v1,&v2,&v3);
        segments->data[2*i] = v1;
        segments->data[2*i+1] = v2;
    }
    
    int j, v,temp;
    for (int i = 0; i<nsegs; i++){
        j = i+1;
        v = segments->data[2*i+1];
        while (j < nsegs){
            if (segments->data[2*j] == v){
                temp = segments->data[2*(i+1)];
                segments->data[2*(i+1)] = segments->data[2*j];
                segments->data[2*j] = temp;

                temp = segments->data[2*(i+1)+1];
                segments->data[2*(i+1)+1] = segments->data[2*j+1];
                segments->data[2*j+1] = temp;
                break;
            }
            j++;
        }

    }

    for(int i = 225; i<nsegs; i++){
        temp = segments->data[2*i];
        segments->data[2*i] = segments->data[2*i+1];
        segments->data[2*i+1] = temp;
    }
    
    fclose(fid);
    return;
}