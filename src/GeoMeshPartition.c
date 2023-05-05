#include "GeoMesh.h"

void stereo_Project_up(struct DoubleMatrix* points2d, struct DoubleMatrix* points3d);
void stereo_Project_down(struct DoubleMatrix* points2d, struct DoubleMatrix* points3d);
void centerpoint(struct DoubleMatrix* points3d, double* C, int csample);
void conmap(struct DoubleMatrix* points3d, double* C, struct DoubleMatrix* points3d_map);
void GeoMesh_partition_kernel(struct Mesh* msh, int part);

void GeoMesh_partition(struct Mesh* msh, int type, int npartitions){
    msh->mprts.type = type;
    msh->mprts.npartitions = npartitions;
    int nelems = msh->nelems;
    int nv = msh->coords.nrows;
    int npoints = (type==1)?nv:nelems;

    // creating mesh crs
    if (msh->hasGraph){
        free(msh->grph.row_idx); free(msh->grph.col_idx);
    }
    //Mesh_Graphinit(msh, type);

    // create initial partition (include all elements or nodes)
    msh->mprts.parts_idx = (int*) malloc((npartitions+1)*sizeof(int));
    msh->mprts.parts_idx[0] = 0; msh->mprts.parts_idx[1] = npoints;
    msh->mprts.parts = (int*) malloc(npoints*sizeof(int));
    for (int ii = 0; ii<npoints; ii++){
        msh->mprts.parts[ii] = ii;
    }

    GeoMesh_partition_kernel(msh, 0);
    return;
}

void GeoMesh_partition_kernel(struct Mesh* msh, int part){
    int type = msh->mprts.type;
    int npoints = msh->mprts.parts_idx[part+1]-msh->mprts.parts_idx[part];
    struct DoubleMatrix* points2d = (struct DoubleMatrix*) malloc(sizeof(struct DoubleMatrix));
    struct DoubleMatrix* points3d = (struct DoubleMatrix*) malloc(sizeof(struct DoubleMatrix));
    *points2d = DoubleMatrix_create(npoints, 2);
    *points3d = DoubleMatrix_create(npoints, 3);

    // create iteration numbers
    int ntries = 30;
    int nlines = (int) pow(((double)ntries /2.0), 2.0/3.0);
    int nouter = (int) ceil(log((double) ntries-nlines+1) / log(20.0));
    int ninner = (ntries - nlines) / nouter;
    int csample = Min(npoints, pow(5,4));

    int ptr;
    if (type == 1){
        for (int ii = 0; ii<npoints; ii++){
            ptr = msh->mprts.parts[msh->mprts.parts_idx[part] + ii];
            points2d->data[2*ii] = msh->coords.data[2*ptr];
            points2d->data[2*ii+1] = msh->coords.data[2*ptr+1];
        }
    } else if (type == 2){
        for (int ii = 0; ii<npoints; ii++){
            ptr = msh->mprts.parts[msh->mprts.parts_idx[part] + ii];
            points2d->data[2*ii] = (msh->coords.data[2*msh->elems.data[3*ptr]] + \
            msh->coords.data[2*msh->elems.data[3*ptr+1]] + \
            msh->coords.data[2*msh->elems.data[3*ptr+2]])/3.0;

            points2d->data[2*ii+1] = (msh->coords.data[2*msh->elems.data[3*ptr]+1] + \
            msh->coords.data[2*msh->elems.data[3*ptr+1]+1] + \
            msh->coords.data[2*msh->elems.data[3*ptr+2]+1])/3.0;
        }
    }

    // Rescaling points to be between -1 and 1
    double c[2] = {0.0,0.0};
    for (int ii = 0; ii<npoints; ii++){
        c[0]+= points2d->data[2*ii];
        c[1]+= points2d->data[2*ii+1];
    }
    c[0] = c[0] / (double)npoints;
    c[1] = c[1] / (double)npoints;

    for (int ii = 0; ii<npoints; ii++){
        points2d->data[2*ii]-= c[0];
        points2d->data[2*ii+1]-= c[1];
    }
    double M = 0.0;
    for (int ii = 0; ii<npoints; ii++){
        M = Max(M, Max(Absolute(points2d->data[2*ii]), Absolute(points2d->data[2*ii+1])));
    }

    for (int ii = 0; ii<npoints; ii++){
        points2d->data[2*ii] /= M;
        points2d->data[2*ii+1] /= M;
    }

    // project onto sphere
    stereo_Project_up(points2d, points3d);
    msh->coords = *points3d;

    struct DoubleMatrix* points3d_map = (struct DoubleMatrix*) malloc(sizeof(struct DoubleMatrix));
    *points3d_map = DoubleMatrix_create(npoints, 3);

    // main loop
    double bestcirclequality = 1.0e308;
    double circlequality;
    double* C = (double*) malloc(3*sizeof(double));
    for (int n = 0; n<nouter; n++){
        centerpoint(points3d, C, csample);

        // find conmap (points3d, C, points3d_map)

        // circle circle,quality (double[3], double) 
        // msh.grph, points3d_map, ninner
        if (circlequality < bestcirclequality){
            bestcirclequality = circlequality;
            // others
        }
    }
    free(C);

    // partition and reorganize data


    free(points3d->data); free(points2d->data); free(points3d_map->data);
    free(points3d); free(points2d); free(points3d_map);
    return;
}

void centerpoint(struct DoubleMatrix* points3d, double* C, int csample){
    int npoints = points3d->nrows;
    double npointsd = (double)npoints;
    double ndim = (double)points3d->ncols;
    int ndimi = points3d->ncols;
    double n = Min((double)csample, npointsd);
    n = (ndim+1.0)*floor((n-1)/(ndim+1)) + 1;

    int sz = (int) ceil(n*(1 + 1/(ndim+1)));
    // allocating points3ds
    struct DoubleMatrix* points3ds = (struct DoubleMatrix*) malloc(sizeof(struct DoubleMatrix));
    points3ds->nrows = sz; points3ds->ncols = (int)ndim;
    points3ds->data = (double*) malloc(sz*ndimi*sizeof(double));

    // allocating matrix to find nullvector
    struct DoubleMatrix* A = (struct DoubleMatrix*) malloc(sizeof(struct DoubleMatrix));
    A->nrows = ndimi+2; A->ncols = ndimi+1;
    A->data = (double*) malloc((ndimi+2)*(ndimi+1)*sizeof(double));
    for (int i = 0; i<ndimi+2; i++){
        A->data[(ndimi+1)*i] = 1.0;
    }

    struct DoubleMatrix* Q = (struct DoubleMatrix*) malloc(sizeof(struct DoubleMatrix));
    struct DoubleMatrix* R = (struct DoubleMatrix*) malloc(sizeof(struct DoubleMatrix));

    // random permutation
    int* sample = (int*) malloc(npoints*sizeof(int));
    int j,temp;
    for (int i = 0; i<npoints; i++){sample[i]=i;}
    for (int i = npoints-1; i>=0; --i){
        j = rand() % (i+1);
        temp = sample[i];
        sample[i] = sample[j];
        sample[j] = temp;
    }
    for (int i = 0; i<n; i++){
        for (int j = 0; j<ndimi; j++){
            points3ds->data[3*i+j] = points3d->data[3*sample[i]+ j];
        }
    }

    int qhead = 0;
    int qtail = n-1;
    int pos[ndimi];
    int kk;
    double sum;
    while (qhead<qtail){
        qtail++;
        // store A
        for (int i = 0; i<ndimi+2; i++){
            for (int j = 0; j<ndimi; j++){
                A->data[(ndimi+1)*i + j+1] = \
                points3ds->data[ndimi*(qhead+i) + j];
            }
        }

        // perform radon
        kk = 0;
        sum = 0.0;
        QR(A,Q,R);
        for (int i = 0; i<ndimi+2; i++){
            if (Q->data[i*(ndimi+2) + ndimi+1] > 0){
                pos[kk] = i;
                kk++;
                sum += Q->data[i*(ndimi+2) + ndimi+1];
            }
        }

        // make new qtail
        for (int j = 0; j<ndimi; j++){
            // nullvec(pos)*A(pos,2:4)/sum(nullvec(pos));
            points3ds->data[ndimi*qtail + j] = 0.0;
            for (int k = 0; k<kk; k++){
                points3ds->data[ndimi*qtail + j] += \
                    Q->data[pos[k]*(ndimi+2) + ndimi+1] * \
                    A->data[(ndimi+1)*pos[k] + j+1]/sum;
            }
        }
        free(Q->data); free(R->data);

        qhead+=ndimi+2;
    }
    
    // store answer
    for(int j = 0; j<ndimi; j++){
        C[j] = points3ds->data[ndimi*qtail + j];
    }
    printf("centerpoint: %f %f %f\n",C[0],C[1],C[2]);

    free(sample);
    free(A->data); free(A);
    free(Q); free(R);
    free(points3ds->data);
    free(points3ds);
}

void conmap(struct DoubleMatrix* points3d, double* C, struct DoubleMatrix* points3d_map){
    int ndim = points3d->ncols;
    int npoints = points3d->nrows;
    struct DoubleMatrix* A = (struct DoubleMatrix*) malloc(sizeof(struct DoubleMatrix));
    struct DoubleMatrix* Q = (struct DoubleMatrix*) malloc(sizeof(struct DoubleMatrix));
    struct DoubleMatrix* R = (struct DoubleMatrix*) malloc(sizeof(struct DoubleMatrix));
    struct DoubleMatrix* points_down = (struct DoubleMatrix*) malloc(sizeof(struct DoubleMatrix));

    A->nrows = ndim; A->ncols = 1;
    A->data = (double*) malloc(ndim*sizeof(double));
    for (int i = 0; i<ndim; i++){
        A->data[ndim-1-i] = C[i];
    }

    points_down->nrows = npoints;
    points_down->ncols = ndim-1;
    points_down->data = (double*) malloc(npoints*(ndim-1)*sizeof(double));

    QR(A,Q,R);
    double alpha = R->data[0];
    alpha = sqrt((1.0+alpha)/(1.0-alpha));

    // reverse Q
    int start = 0; int end = (ndim)*(ndim)-1; int temp;
    while (start < end){
        temp = Q->data[start];
        Q->data[start] = Q->data[end];
        Q->data[end] = temp;
        start++; end--;
    }

    // points3d_map = points3d*Q
    for (int i = 0; i<npoints; i++){
        for (int j = 0; j<ndim; j++){
            points3d_map->data[ndim*i + j] = 0.0;
            for (int k = 0; k<ndim; k++){
                points3d_map->data[ndim*i + j] += \
                    points3d->data[ndim*i + k] * \
                    Q->data[ndim*k + j];
            }
        }
    }

    // project down

    free(points_down->data); free(points_down);
    free(A->data); free(A);
    free(Q->data); free(Q);
    free(R->data); free(R);
    return;
}




void stereo_Project_up(struct DoubleMatrix* points2d, struct DoubleMatrix* points3d){
    int npoints = points2d->nrows;
    assert(points3d->nrows == npoints);
    double nrmsqp1;
    for (int ii = 0; ii<npoints; ii++){
        nrmsqp1 = 1.0 + points2d->data[2*ii]*points2d->data[2*ii] + \
        points2d->data[2*ii+1]*points2d->data[2*ii+1];
        points3d->data[3*ii] = 2.0*points2d->data[2*ii] / nrmsqp1;
        points3d->data[3*ii+1] = 2.0*points2d->data[2*ii+1] / nrmsqp1;
        points3d->data[3*ii+2] = 1.0 - 2.0 / nrmsqp1;
    }
    return;
}

void stereo_Project_down(struct DoubleMatrix* points2d, struct DoubleMatrix* points3d){
    int npoints = points2d->nrows;
    
}