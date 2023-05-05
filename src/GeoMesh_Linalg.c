#include "GeoMesh.h"

void QR(const struct DoubleMatrix* A, struct DoubleMatrix* Q, struct DoubleMatrix* R){
    int m = A->nrows;
    int n = A->ncols;
    *Q = Eye(m);
    *R = DoubleMatrix_create(m,n);
    for (int i = 0; i<m*n; i++){
        R->data[i] = A->data[i];
    }

    struct DoubleMatrix* Rbuff = (struct DoubleMatrix*) malloc(sizeof(struct DoubleMatrix));
    Rbuff->nrows = m; Rbuff->ncols = n;
    Rbuff->data = (double*) malloc(m*n*sizeof(double));
    struct DoubleMatrix* Qbuff = (struct DoubleMatrix*) malloc(sizeof(struct DoubleMatrix));
    Qbuff->nrows = m; Qbuff->ncols = m;
    Qbuff->data = (double*) malloc(m*m*sizeof(double));

    int sz = (m>=n)?m:n;
    double* r = (double*) malloc(sz*sizeof(double));
    double* w = (double*) malloc(sz*sizeof(double));
    int i,k;
    double normr, s, tau, in, u1;
    for (int j = 0; j<n; j++){ // loop through columns
        sz = m-j;

        // R matrix buffer
        for (k = 0; k<sz; k++){
            r[k] = R->data[(k+j)*n + j];
            for (i = 0; i<n; i++){
                Rbuff->data[k*n+i] = R->data[(k+j)*n + i];
            }
        }

        // Q matrix buffer
        for (k = 0; k<m; k++){
            for (i = 0; i<sz; i++){
                Qbuff->data[k*m + i] = Q->data[k*m + (i+j)];
            }
        }

        normr = 0.0;
        for (k=0; k<sz;k++){normr+=r[k]*r[k];}
        normr = sqrt(normr);
        s = Sign(R->data[n*j + j]);
        u1 = R->data[j*n+j] + s*normr;
        w[0] = 1.0; for (k=1; k<sz; k++){w[k]=r[k]/u1;}
        tau = s*u1/normr;
        in = 0.0;
        for (k=0; k<sz; k++){in += r[k]*w[k];}
        for (k=0; k<sz; k++){r[k] -= tau*in*w[k];}

        // use r to store w*Rbuff
        for (k = 0; k<n; k++){
            r[k] = 0.0;
            for (i=0; i<sz; i++){
                r[k] += w[i]*Rbuff->data[i*n + k];
            }
        }

        // updating R
        for (k = 0; k<sz; k++){
            for (i = 0; i<n; i++){
                R->data[(k+j)*n + i] -= tau*w[k]*r[i];
            }
        }

        // use r to store Qbuff*w
        for (k = 0; k<m; k++){
            r[k] = 0.0;
            for (i=0; i<sz; i++){
                r[k] += Qbuff->data[k*m + i]*w[i];
            }
        }

        // updating Q
        for (k=0; k<m; k++){
            for (i=0; i<sz; i++){
                Q->data[k*m + (i+j)] -= r[k]*tau*w[i];
            }
        }
    }
    
    free(Qbuff->data); free(Qbuff);
    free(Rbuff->data); free(Rbuff);
    free(r); free(w);
    return;
}

struct DoubleMatrix Eye(int n){
    struct DoubleMatrix I = Zeros(n,n);
    for (int i = 0; i<n; i++){
        I.data[I.ncols*i + i] = 1;
    }
    return I;
}

struct DoubleMatrix Zeros(int m, int n){
    struct DoubleMatrix A = DoubleMatrix_create(m,n);
    memset(A.data, 0.0, m * n * sizeof(double));
    return A;
}