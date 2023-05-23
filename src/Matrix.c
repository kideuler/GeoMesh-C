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

void Mesh_resize(struct Mesh* msh, int newsz){
    DoubleMatrix_resize(&msh->coords,newsz);
    IntMatrix_resize(&msh->elems,newsz);
    IntMatrix_resize(&msh->sibhfs,newsz);
    
    bool* delete_elem = (bool*) malloc(newsz*sizeof(bool));
    memset(delete_elem, false, newsz*sizeof(bool));
    for(int i = 0; i<msh->nelems; i++){
        delete_elem[i] = msh->delete_elem[i];
    }
    free(msh->delete_elem);
    msh->delete_elem = delete_elem;
    
    bool* on_bdy = (bool*) malloc(newsz*sizeof(bool));
    memset(on_bdy, false, newsz*sizeof(bool));
    for(int i = 0; i<msh->nelems; i++){
        on_bdy[i] = msh->on_boundary[i];
    }
    free(msh->on_boundary);
    msh->on_boundary = on_bdy;
    return;
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

struct DoubleMatrix prod(struct DoubleMatrix* A, struct DoubleMatrix* B){
    int m = A->nrows;
    int r = A->ncols;
    int n = B->ncols;

    struct DoubleMatrix C = DoubleMatrix_create(m,n);
    for (int i = 0; i<m; i++){
        for (int j = 0; j<n; j++){
            C.data[n*i+j] = 0.0;
            for (int k = 0; k<r; k++){
                C.data[n*i+j] += A->data[r*i + k]*B->data[n*k + j];
            }
        }
    }
    return C;
}
struct DoubleMatrix subtract(struct DoubleMatrix* A, struct DoubleMatrix* B){
    int m = A->nrows;
    int n = B->ncols;

    struct DoubleMatrix C = DoubleMatrix_create(m,n);
    for (int i = 0; i<m; i++){
        for (int j = 0; j<n; j++){
            C.data[n*i+j] = A->data[n*i+j] - B->data[n*i+j];
        }
    }
    return C;
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

void load_monalisa(struct DoubleMatrix* coords){
    FILE *fid;
    int nv,nsegs,v1,v2,v3;
    float x, y, z;
    free(coords->data);

    fid = fopen("data/monalisa.dat","r");

    fscanf(fid,"%d",&nv);
    *coords = DoubleMatrix_create(nv,2);

    for(int i=0; i<nv; i++){
        fscanf(fid,"%d %f %f", &v1,&x,&y);
        coords->data[2*i] = x;
        coords->data[2*i+1] = y;
    }

    fclose(fid);
    return;
}

double Circle(struct IntMatrix* segments, struct DoubleMatrix* coords, int npoints, bool box){
    int size  = box?(npoints+4):npoints;
    *segments = IntMatrix_create(npoints,2);
    *coords = DoubleMatrix_create(size,2);

    double t;
    for (int i = 0; i<npoints; i++){
        t = 2*M_PI*(((double) i) / ((double) npoints));
        coords->data[2*i] = 0.5*cos(t)+0.5;
        coords->data[2*i+1] = 0.5*sin(t)+0.5;
    }

    if (box){
        coords->data[2*npoints] = 0.0001;
        coords->data[2*npoints+1] = 0.0;
        coords->data[2*(npoints+1)] = 1.0;
        coords->data[2*(npoints+1)+1] = 0.0;
        coords->data[2*(npoints+2)] = 1.0;
        coords->data[2*(npoints+2)+1] = 1.0;
        coords->data[2*(npoints+3)] = 0.0;
        coords->data[2*(npoints+3)+1] = 1.0;
    }

    double h = 0.0;
    if (box) {
        for(int i = 0; i<npoints; i++){
            segments->data[2*i+1]=i; segments->data[2*i]=(i+1)%npoints;
            h += sqrt(pow(coords->data[2*((i+1)%npoints)]-coords->data[2*i],2) + \
            pow(coords->data[2*((i+1)%npoints)+1]-coords->data[2*i+1],2));
        }
    } else {
        for(int i = 0; i<npoints; i++){
            segments->data[2*i]=i; segments->data[2*i+1]=(i+1)%npoints;
            h += sqrt(pow(coords->data[2*((i+1)%npoints)]-coords->data[2*i],2) + \
            pow(coords->data[2*((i+1)%npoints)+1]-coords->data[2*i+1],2));
        }
    }

    h = h/((double)npoints);
    return h;
}

double Ellipse(struct IntMatrix* segments, struct DoubleMatrix* coords, int npoints, bool box){
    int size  = box?(npoints+4):npoints;
    *segments = IntMatrix_create(npoints,2);
    *coords = DoubleMatrix_create(size,2);

    double t;
    for (int i = 0; i<npoints; i++){
        t = 2*M_PI*(((double) i) / ((double) npoints));
        coords->data[2*i] = 0.2*cos(t)+0.5;
        coords->data[2*i+1] = 0.1*sin(t)+0.5;
    }

    if (box){
        coords->data[2*npoints] = 0.0001;
        coords->data[2*npoints+1] = 0.0;
        coords->data[2*(npoints+1)] = 1.0;
        coords->data[2*(npoints+1)+1] = 0.0;
        coords->data[2*(npoints+2)] = 1.0;
        coords->data[2*(npoints+2)+1] = 1.0;
        coords->data[2*(npoints+3)] = 0.0;
        coords->data[2*(npoints+3)+1] = 1.0;
    }

    double h = 0.0;
    if (box) {
        for(int i = 0; i<npoints; i++){
            segments->data[2*i+1]=i; segments->data[2*i]=(i+1)%npoints;
            h += sqrt(pow(coords->data[2*((i+1)%npoints)]-coords->data[2*i],2) + \
            pow(coords->data[2*((i+1)%npoints)+1]-coords->data[2*i+1],2));
        }
    } else {
        for(int i = 0; i<npoints; i++){
            segments->data[2*i]=i; segments->data[2*i+1]=(i+1)%npoints;
            h += sqrt(pow(coords->data[2*((i+1)%npoints)]-coords->data[2*i],2) + \
            pow(coords->data[2*((i+1)%npoints)+1]-coords->data[2*i+1],2));
        }
    }

    h = h/((double)npoints);
    return h;
}

double Flower(struct IntMatrix* segments, struct DoubleMatrix* coords, int npoints, bool box){
    int size  = box?(npoints+4):npoints;
    *segments = IntMatrix_create(npoints,2);
    *coords = DoubleMatrix_create(size,2);

    double t;
    for (int i = 0; i<npoints; i++){
        t = 2*M_PI*(((double) i) / ((double) npoints));
        coords->data[2*i] = (0.25 + 0.1*sin(5*t))*cos(t)+0.5;
        coords->data[2*i+1] = (0.25 + 0.1*sin(5*t))*sin(t)+0.5;
    }

    if (box){
        coords->data[2*npoints] = 0.0001;
        coords->data[2*npoints+1] = 0.0;
        coords->data[2*(npoints+1)] = 1.0;
        coords->data[2*(npoints+1)+1] = 0.0;
        coords->data[2*(npoints+2)] = 1.0;
        coords->data[2*(npoints+2)+1] = 1.0;
        coords->data[2*(npoints+3)] = 0.0;
        coords->data[2*(npoints+3)+1] = 1.0;
    }

    double h = 0.0;
    if (box) {
        for(int i = 0; i<npoints; i++){
            segments->data[2*i+1]=i; segments->data[2*i]=(i+1)%npoints;
            h += sqrt(pow(coords->data[2*((i+1)%npoints)]-coords->data[2*i],2) + \
            pow(coords->data[2*((i+1)%npoints)+1]-coords->data[2*i+1],2));
        }
    } else {
        for(int i = 0; i<npoints; i++){
            segments->data[2*i]=i; segments->data[2*i+1]=(i+1)%npoints;
            h += sqrt(pow(coords->data[2*((i+1)%npoints)]-coords->data[2*i],2) + \
            pow(coords->data[2*((i+1)%npoints)+1]-coords->data[2*i+1],2));
        }
    }

    h = h/((double)npoints);
    return h;
}

double Airfoil(struct IntMatrix* segments, struct DoubleMatrix* coords, int npoints, bool box){
    int size = 2*npoints;
    size  = box?(size+4):size;
    *segments = IntMatrix_create(2*npoints,2);
    *coords = DoubleMatrix_create(size,2);
    double t;
    double yc,theta, yt;
    double p = 0.4;
    double m = 0.02;
    double T = 0.12;
    int nv = 0;
    for (int i = 0; i<=npoints; i++){
        t = (((double) i) / ((double) npoints));
        yt = 5*T*(0.2969*sqrt(t) - 0.1260*t - 0.3516*t*t + 0.2863*t*t*t - 0.1015*t*t*t*t);
        yc = (t<=p) ? (m/pow(p,2))*(2*p*t - t*t) : (m/pow(1-p,2))*((1-2*p)+2*p*t-t*t);
        theta = atan((t<=p) ? 2*(m/pow(p,2))*(p-t) : 2*(m/pow(1-p,2))*(p-t));
        coords->data[2*nv] = 0.5*(t - yt*sin(theta))+0.25;
        coords->data[2*nv+1] = 0.5*(yc + yt*cos(theta))+0.5;
        nv++;
    }
    for (int i = 1; i<npoints; i++){
        t = 1-(((double) i) / ((double) npoints));
        yt = 5*T*(0.2969*sqrt(t) - 0.1260*t - 0.3516*t*t + 0.2863*t*t*t - 0.1015*t*t*t*t);
        yc = (t<=p) ? (m/pow(p,2))*(2*p*t - t*t) : (m/pow(1-p,2))*((1-2*p)+2*p*t-t*t);
        theta = atan((t<=p) ? 2*(m/pow(p,2))*(p-t) : 2*(m/pow(1-p,2))*(p-t));
        coords->data[2*(nv)] = 0.5*(t + yt*sin(theta))+0.25;
        coords->data[2*(nv)+1] = 0.5*(yc - yt*cos(theta))+0.5;
        nv++;
    }


    double h = 0.0;
    npoints = 2*npoints;

    if (box){
        coords->data[2*npoints] = 0.0001;
        coords->data[2*npoints+1] = 0.0;
        coords->data[2*(npoints+1)] = 1.0;
        coords->data[2*(npoints+1)+1] = 0.0;
        coords->data[2*(npoints+2)] = 1.0;
        coords->data[2*(npoints+2)+1] = 1.0;
        coords->data[2*(npoints+3)] = 0.0;
        coords->data[2*(npoints+3)+1] = 1.0;
    }

    if (!box) {
        for(int i = 0; i<npoints; i++){
            segments->data[2*i+1]=i; segments->data[2*i]=(i+1)%npoints;
            h += sqrt(pow(coords->data[2*((i+1)%npoints)]-coords->data[2*i],2) + \
            pow(coords->data[2*((i+1)%npoints)+1]-coords->data[2*i+1],2));
        }
    } else {
        for(int i = 0; i<npoints; i++){
            segments->data[2*i]=i; segments->data[2*i+1]=(i+1)%npoints;
            h += sqrt(pow(coords->data[2*((i+1)%npoints)]-coords->data[2*i],2) + \
            pow(coords->data[2*((i+1)%npoints)+1]-coords->data[2*i+1],2));
        }
    }

    h = h/((double)npoints);
    return h;
}


double Cube(struct DoubleMatrix* coords, int npoints){
    *coords = DoubleMatrix_create(npoints*npoints*npoints,3);

    int ii,jj,kk,nn;
    nn=0;
    double h = 1.0 / ((double)npoints-1.0);
    for (ii=0;ii<npoints;ii++){
        for (jj=0;jj<npoints;jj++){
            for (kk=0;kk<npoints;kk++){
                coords->data[3*nn] = h*((double) kk);
                coords->data[3*nn+1] = h*((double) jj);
                coords->data[3*nn+2] = h*((double) ii);
                nn++;
            }
        }
    }
}

double Sphere(struct DoubleMatrix* coords, int npoints, bool box){
    int size  = box?(npoints+8):npoints;
    *coords = DoubleMatrix_create_Random(size,3,-M_PI,M_PI);

    double phi = M_PI*(sqrt(5.0)-1);
    double x,y,z,radius,theta;
    for (int i = 0; i<size; i++){
        y = 1.0 - (((double) i) / ((double) size-1))*2.0;
        radius = sqrt(1.0-y*y);
        theta = phi * ((double) i);
        x = cos(theta)*radius;
        z = sin(theta)*radius;

        coords->data[3*i] = 0.3*x+0.5;
        coords->data[3*i+1] = 0.3*y+0.5;
        coords->data[3*i+2] = 0.3*z+0.5;
    }

    if (box){
        coords->data[3*(npoints)] = 0.0;
        coords->data[3*(npoints)+1] = 0.0;
        coords->data[3*(npoints)+2] = 0.0;

        coords->data[3*(npoints+1)] = 1.0;
        coords->data[3*(npoints+1)+1] = 0.0;
        coords->data[3*(npoints+1)+2] = 0.0;

        coords->data[3*(npoints+2)] = 1.0;
        coords->data[3*(npoints+2)+1] = 1.0;
        coords->data[3*(npoints+2)+2] = 0.0;

        coords->data[3*(npoints+3)] = 0.0;
        coords->data[3*(npoints+3)+1] = 1.0;
        coords->data[3*(npoints+3)+2] = 0.0;

        coords->data[3*(npoints+4)] = 0.0;
        coords->data[3*(npoints+4)+1] = 0.0;
        coords->data[3*(npoints+4)+2] = 1.0;

        coords->data[3*(npoints+5)] = 1.0;
        coords->data[3*(npoints+5)+1] = 0.0;
        coords->data[3*(npoints+5)+2] = 1.0;

        coords->data[3*(npoints+6)] = 1.0;
        coords->data[3*(npoints+6)+1] = 1.0;
        coords->data[3*(npoints+6)+2] = 1.0;

        coords->data[3*(npoints+7)] = 0.0;
        coords->data[3*(npoints+7)+1] = 1.0;
        coords->data[3*(npoints+7)+2] = 1.0;

    }

    double h = 0.0;
    return h;
}

double Ellipsoid(struct DoubleMatrix* coords, int npoints, bool box){
    int size  = box?(npoints+8):npoints;
    *coords = DoubleMatrix_create_Random(size,3,-M_PI,M_PI);

    double phi = M_PI*(sqrt(5.0)-1);
    double x,y,z,radius,theta;
    for (int i = 0; i<size; i++){
        y = 1.0 - (((double) i) / ((double) size-1))*2.0;
        radius = sqrt(1.0-y*y);
        theta = phi * ((double) i);
        x = cos(theta)*radius;
        z = sin(theta)*radius;

        coords->data[3*i] = (7.0/20.0)*x+0.5;
        coords->data[3*i+1] = 0.25*y+0.5;
        coords->data[3*i+2] = 0.2*z+0.5;
    }

    if (box){
        coords->data[3*(npoints)] = 0.0;
        coords->data[3*(npoints)+1] = 0.0;
        coords->data[3*(npoints)+2] = 0.0;

        coords->data[3*(npoints+1)] = 1.0;
        coords->data[3*(npoints+1)+1] = 0.0;
        coords->data[3*(npoints+1)+2] = 0.0;

        coords->data[3*(npoints+2)] = 1.0;
        coords->data[3*(npoints+2)+1] = 1.0;
        coords->data[3*(npoints+2)+2] = 0.0;

        coords->data[3*(npoints+3)] = 0.0;
        coords->data[3*(npoints+3)+1] = 1.0;
        coords->data[3*(npoints+3)+2] = 0.0;

        coords->data[3*(npoints+4)] = 0.0;
        coords->data[3*(npoints+4)+1] = 0.0;
        coords->data[3*(npoints+4)+2] = 1.0;

        coords->data[3*(npoints+5)] = 1.0;
        coords->data[3*(npoints+5)+1] = 0.0;
        coords->data[3*(npoints+5)+2] = 1.0;

        coords->data[3*(npoints+6)] = 1.0;
        coords->data[3*(npoints+6)+1] = 1.0;
        coords->data[3*(npoints+6)+2] = 1.0;

        coords->data[3*(npoints+7)] = 0.0;
        coords->data[3*(npoints+7)+1] = 1.0;
        coords->data[3*(npoints+7)+2] = 1.0;

    }

    double h = 0.0;
    return h;
}