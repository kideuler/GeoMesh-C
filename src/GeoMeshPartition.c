#include "GeoMesh.h"

void stereo_Project_up(struct DoubleMatrix* points2d, struct DoubleMatrix* points3d);
void stereo_Project_down(struct DoubleMatrix* points2d, struct DoubleMatrix* points3d);
void partition(int* part1, int* npart1, int* part2, int* npart2, double* quality, int qual, double percent_cut);
void centerpoint(struct DoubleMatrix* points3d, double* C, int csample);
void conmap(struct DoubleMatrix* points3d, double* C, struct DoubleMatrix* points3d_map);
int separate_circle(struct Mesh* msh, struct DoubleMatrix* points_map, int ntries, double* Circle, int part, double percent_cut);
int separate_line(struct Mesh* msh, struct DoubleMatrix* points_map, int ntries, double* line, int part, double percent_cut);
void GeoMesh_partition_kernel(struct Mesh* msh, int part, double percent_cut);
void Geomesh_partitionbase2(struct Mesh* msh, int part, int nparts);

// quicksort private functions
void quicksort_kernel(double *x, int lo, int hi);
int quicksort_partition(double *x, int lo, int hi);

void GeoMesh_partition(struct Mesh* msh, int type, int npartitions){
    msh->mprts.type = type;
    msh->mprts.npartitions = 1;
    msh->hasPartition = true;
    int nelems = msh->nelems;
    int nv = msh->coords.nrows;
    int npoints = (type==1)?nv:nelems;

    // creating mesh crs
    if (msh->hasGraph){
        //free(msh->grph.row_idx); free(msh->grph.col_idx);
    }
    //Mesh_Graphinit(msh, type);

    // create initial partition (include all elements or nodes)
    msh->mprts.parts_idx = (int*) malloc((npoints)*sizeof(int));
    msh->mprts.parts_idx[0] = 0; msh->mprts.parts_idx[1] = npoints;
    msh->mprts.parts = (int*) malloc(npoints*sizeof(int));
    for (int ii = 0; ii<npoints; ii++){
        msh->mprts.parts[ii] = ii;
    }

    if (floor(log2((double)npartitions)) == ceil(log2((double)npartitions))){
        Geomesh_partitionbase2(msh, 0, npartitions);
    } else {
        // repeatedly subtract off portions of highest logs
        int remainder = npartitions;
        int curr;
        double percent_cut;
        while (remainder > 0){
            if (remainder == 2){Geomesh_partitionbase2(msh, 0, 2); break;}
            curr = (int) pow(2.0,floor(log2((double)remainder)));
            percent_cut = 1.0 - (double)curr / (double)remainder;
            GeoMesh_partition_kernel(msh, 0, percent_cut);
            Geomesh_partitionbase2(msh, 1, curr);
            remainder -= curr;
            if (remainder <= 1){break;}
        }
    }
    return;
}

void Geomesh_partitionbase2(struct Mesh* msh, int part, int nparts){
    int niters = (int) log2((double) nparts);
    assert(niters == (int)ceil(log2((double)nparts)));

    for (int n = 0; n<niters; n++){
        for (int i = 0; i<(int)pow(2.0,(double)n+1); i=i+2){
            GeoMesh_partition_kernel(msh, part+i,0.5);
        }
    }
    return;
}

void quicksort(double* x, int sz){
    quicksort_kernel(x,0,sz-1);
}

void GeoMesh_partition_kernel(struct Mesh* msh, int part, double percent_cut){
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

    struct DoubleMatrix* points3d_map = (struct DoubleMatrix*) malloc(sizeof(struct DoubleMatrix));
    *points3d_map = DoubleMatrix_create(npoints, 3);

    // main loop for circlequality
    int bestcirclequality = INT32_MAX;
    int circlequality;
    double* C = (double*) malloc(3*sizeof(double));
    double* BestC = (double*) malloc(3*sizeof(double));
    double* Circle = (double*) malloc(3*sizeof(double));
    double* BestCircle = (double*) malloc(3*sizeof(double));
    for (int n = 0; n<nouter; n++){
        // find centerpoint C
        centerpoint(points3d, C, csample);
        //printf("centerpoint: %f %f %f\n",C[0],C[1],C[2]);

        // map 3d points to centerpoint
        conmap(points3d, C, points3d_map);

        circlequality = separate_circle(msh, points3d_map, ninner, Circle, part, percent_cut);
        if (circlequality < bestcirclequality){
            bestcirclequality = circlequality;
            for (int j = 0; j<3; j++){BestCircle[j] = Circle[j];}
            for (int j = 0; j<3; j++){BestC[j] = C[j];}
        }
    }
    free(C);
    free(Circle);

    // partition by line
    double* L = (double*) malloc(2*sizeof(double));
    int Linequality = separate_line(msh, points2d, nlines, L, part, percent_cut);
    
    // partition and reorganize data
    double* Dotprod = (double*) malloc(npoints*sizeof(double));
    memset(Dotprod, 0.0, npoints*sizeof(double));

    if (bestcirclequality < Linequality) {
        conmap(points3d, BestC, points3d_map);
        for (int i=0;i<npoints;i++){
            for (int j = 0; j<3; j++){
                Dotprod[i] += points3d_map->data[3*i + j]*BestCircle[j];
            }
        }
    } else {
        for (int i=0;i<npoints;i++){
            for (int j = 0; j<2; j++){
                Dotprod[i] += points2d->data[2*i + j]*L[j];
            }
        }
    }

    int* partsnew = (int*) malloc(npoints*sizeof(int));
    int* part1 = (int*) malloc(npoints*sizeof(int)); int npart1;
    int* part2 = (int*) malloc(npoints*sizeof(int)); int npart2;
    partition(part1, &npart1, part2, &npart2, Dotprod, npoints, percent_cut);

    int kk = 0;
    for (int i = 0; i<npart1; i++){
        ptr = msh->mprts.parts[msh->mprts.parts_idx[part] + part1[i]];
        partsnew[kk] = ptr;
        kk++;
    }
    for (int i = 0; i<npart2; i++){
        ptr = msh->mprts.parts[msh->mprts.parts_idx[part] + part2[i]];
        partsnew[kk] = ptr;
        kk++;
    }

    for(int i = 0; i<npoints; i++){
        msh->mprts.parts[msh->mprts.parts_idx[part] + i] = partsnew[i];
    }

    for (int i = msh->mprts.npartitions+1; i>part; i--){
        msh->mprts.parts_idx[i] = msh->mprts.parts_idx[i-1];
    }
    msh->mprts.parts_idx[part+1] = npart1+msh->mprts.parts_idx[part]; 
    msh->mprts.npartitions++;

    free(part1); free(part2); free(L);
    free(partsnew);
    free(Dotprod);
    free(BestCircle); free(BestC);
    free(points3d->data);
    free(points2d->data); 
    free(points3d_map->data);
    free(points3d);
    free(points2d);
    free(points3d_map);
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
    int start = 0; int end = (ndim)*(ndim)-1; double temp;
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
    stereo_Project_down(points_down, points3d_map);

    for (int i = 0; i<npoints; i++){
        for (int j = 0; j<ndim-1;j++){
            points_down->data[(ndim-1)*i + j] = points_down->data[(ndim-1)*i + j]/alpha;
        }
    }

    stereo_Project_up(points_down, points3d_map);

    free(points_down->data); free(points_down);
    free(A->data); free(A);
    free(Q->data); free(Q);
    free(R->data); free(R);
    return;
}

int separate_circle(struct Mesh* msh, struct DoubleMatrix* points_map, int ntries, double* Circle, int part, double percent_cut){
    int npoints = points_map->nrows;
    int ndim = points_map->ncols;

    double* P = (double*) malloc(ndim*ndim*sizeof(double));
    memset(P,0.0,ndim*ndim*sizeof(double));
    for (int i = 0; i<ndim; i++){
        for (int j = 0; j<ndim; j++){
            for (int k = 0; k<npoints; k++){
                P[ndim*i + j] += points_map->data[ndim*k + i]*points_map->data[ndim*k + j];
            }
        }
    }

    double* M = (double*) malloc(ndim*ndim*sizeof(double));
    memset(M,0.0,ndim*ndim*sizeof(double));
    for (int i = 0; i<ndim; i++){
        for (int j = 0; j<ndim; j++){
            for (int k = 0; k<ndim; k++){
                M[ndim*i + j] += P[ndim*i + k]*P[ndim*k + j];
            }
        }
    }

    // setting up quality values
    double* qual = (double*) malloc(npoints*sizeof(double));
    int quality = INT32_MAX;
    int q;
    double rvec[ndim], rcircle[ndim]; 
    double nrm;

    // setting up partition arrays
    int* part1 = (int*) malloc(npoints*sizeof(int)); int npart1;
    int* part2 = (int*) malloc(npoints*sizeof(int)); int npart2;

    for (int n = 0; n<ntries; n++){

        // creating random vector
        for (int j=0;j<ndim;j++){rvec[j]=drand(-1.0,1.0);}

        // converting to random circle
        for (int i=0;i<ndim;i++){rcircle[i]=0.0; for (int j=0;j<ndim;j++){rcircle[i]+=rvec[j]*M[ndim*j+i];}}

        nrm=0.0;
        for (int i = 0; i<ndim; i++){nrm+=rcircle[i]*rcircle[i];}
        nrm = sqrt(nrm);

        // check that norm is not zero
        if (nrm > 1e-8){

            // fill in quality array
            for (int i = 0; i<npoints; i++){qual[i] = 0.0;
                for (int j = 0; j<ndim; j++){
                    qual[i] += points_map->data[ndim*i + j]*rcircle[j];
                }
            }

            // partition the quality
            partition(part1, &npart1, part2, &npart2, qual, npoints, percent_cut);

            // find the cutsize from nnz of modified crs array
            q = MeshGraph_nnz(msh, part1, npart1, part2, npart2, part);
            if (q<quality){
                quality = q;
                for (int j =0; j<ndim; j++){
                    Circle[j] = rcircle[j];
                }
            }
        }
    }


    free(P); free(M); free(qual);
    free(part1); free(part2);

    return quality;
}

int separate_line(struct Mesh* msh, struct DoubleMatrix* points_map, int ntries, double* line, int part, double percent_cut){
    int npoints = points_map->nrows;
    int ndim = 2;

    double* P = (double*) malloc(ndim*ndim*sizeof(double));
    memset(P,0.0,ndim*ndim*sizeof(double));
    for (int i = 0; i<ndim; i++){
        for (int j = 0; j<ndim; j++){
            for (int k = 0; k<npoints; k++){
                P[ndim*i + j] += points_map->data[ndim*k + i]*points_map->data[ndim*k + j];
            }
        }
    }

    double T = P[0]+P[3];
    double D = P[0]*P[3] - P[1]*P[2];
    double eig[2] = {T/2 + sqrt(T*T/4.0 - D), T/2 - sqrt(T*T/4.0 - D)};

    double V[4];
    double nrm1, nrm2;
    if (Absolute(P[2]) > 1e-10){
        nrm1 = sqrt(pow(eig[0]-P[3],2.0) + pow(P[2],2.0));
        nrm2 = sqrt(pow(eig[1]-P[3],2.0) + pow(P[2],2.0));
        V[0] = (eig[0] - P[3])/nrm1;
        V[1] = (eig[1] - P[3])/nrm2;
        V[2] = P[2]/nrm1;
        V[3] = P[2]/nrm2;
    } else if (Absolute(P[1]) > 1e-10){
        nrm1 = sqrt(pow(eig[0]-P[0],2.0) + pow(P[1],2.0));
        nrm2 = sqrt(pow(eig[1]-P[0],2.0) + pow(P[1],2.0));
        V[0] = P[1]/nrm1;
        V[1] = P[1]/nrm2;
        V[2] = (eig[0] - P[0])/nrm1;
        V[3] = (eig[1] - P[0])/nrm2;
    } else {
        V[0] = 1.0; V[1] = 0.0; V[2] = 0.0; V[3] = 1.0;
    }
    double xp = 2*((double)ndim+1)/((double)ntries-1);
    eig[0] = pow(sqrt(eig[0]),xp); eig[1] = pow(sqrt(eig[1]),xp);

    double W[4];
    W[0] = eig[0]*V[0]*V[0] + eig[1]*V[1]*V[1];
    W[1] = eig[0]*V[0]*V[2] + eig[1]*V[1]*V[3];
    W[2] = eig[0]*V[0]*V[2] + eig[1]*V[1]*V[3];
    W[3] = eig[0]*V[2]*V[2] + eig[1]*V[3]*V[3];

    double* qual = (double*) malloc(npoints*sizeof(double));
    int quality = INT32_MAX;
    int* part1 = (int*) malloc(npoints*sizeof(int)); int npart1;
    int* part2 = (int*) malloc(npoints*sizeof(int)); int npart2;

    double rvec[2];
    double rline[2];
    double nrml;
    int q;
    for (int n = 0; n<ntries; n++){
        rvec[0] = drand(-1.0,1.0);
        rvec[1] = drand(-1.0,1.0);

        rline[0] = rvec[0]*W[0] + rvec[1]*W[2];
        rline[1] = rvec[0]*W[1] + rvec[1]*W[3];
        nrml = sqrt(rline[0]*rline[0] + rline[1]*rline[1]);
        rline[0] = rline[0]/nrml;
        rline[1] = rline[1]/nrml;
        

        for (int i = 0; i<npoints; i++){qual[i] = 0.0;
            for (int j = 0; j<ndim; j++){
                qual[i] += points_map->data[ndim*i + j]*rline[j];
            }
        }
        partition(part1, &npart1, part2, &npart2, qual, npoints, percent_cut);

        q = MeshGraph_nnz(msh, part1, npart1, part2, npart2, part);
        if (q<quality){
            quality = q;
            for (int j =0; j<ndim; j++){
                line[j] = rline[j];
            }
        }
    }
    
    free(P); free(qual);
    free(part1); free(part2);
}

void partition(int* part1, int* npart1, int* part2, int* npart2, double* quality, int nqual, double percent_cut){
    *npart1 = 0; *npart2 = 0;
    // create buffer for sorting
    double* sortbuffer = (double*) malloc(nqual*sizeof(double));
    for(int n=0;n<nqual;n++){sortbuffer[n]=quality[n];}

    // sort quality
    quicksort(sortbuffer, nqual);
    double split;

    // find upper and lower 
    int upper = (int) ceil(((double)(nqual-1))*percent_cut);
    int lower = (int) floor(((double)(nqual-1))*percent_cut);
    if (upper == lower){split = sortbuffer[upper];}
    else{split = (sortbuffer[lower]+sortbuffer[upper])/2.0;}
    int kk = 0;
    for (int n = 0; n<nqual; n++){
        if (quality[n] < split){
            part1[*npart1] = n;
            *npart1 = *npart1 + 1;
        } else if (quality[n] > split) {
            part2[*npart2] = n;
            *npart2 = *npart2 + 1;
        } else {
            if (kk%2 == 0){
                part1[*npart1] = n;
                *npart1 = *npart1 + 1;
            }else{
                part2[*npart2] = n;
                *npart2 = *npart2 + 1;
            }
            kk++;
        }
    }
    free(sortbuffer);
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
    int npoints = points3d->nrows;
    int ndim = points3d->ncols;
    for (int i = 0; i<npoints; i++){
        for (int j = 0; j<ndim-1; j++){
            points2d->data[(ndim-1)*i + j] = points3d->data[ndim*i + j] / \
                (1.0 - points3d->data[ndim*i + (ndim-1)]);
        }
    }
    return;
}

void quicksort_kernel(double *x, int lo, int hi){
    if (lo >= hi || lo < 0){
        return;
    }

    int p = quicksort_partition(x,lo,hi);
    quicksort_kernel(x,lo,p-1);
    quicksort_kernel(x,p+1,hi);
    return;
}

int quicksort_partition(double *x, int lo, int hi){
    double pivot = x[hi];

    int i = lo - 1;
    double temp;
    for (int j = lo; j<hi; j++){
        if (x[j] <= pivot){
            i = i+1;
            temp = x[i];
            x[i] = x[j];
            x[j] = temp;
        }
    }

    i=i+1;
    temp = x[i];
    x[i] = x[hi];
    x[hi] = temp;
    return i;
}

int MeshGraph_nnz(struct Mesh* msh, int* a, int sza, int* b, int szb, int part){
    int npoints = (msh->mprts.type == 1)?msh->coords.nrows:msh->nelems;
    int j,ptr1,ptr2;
    bool* mask = (bool*) malloc(npoints*sizeof(bool));
    memset(mask, false, npoints*sizeof(bool));
    for (int n = 0; n<szb; n++){
        ptr2 = msh->mprts.parts[msh->mprts.parts_idx[part]+b[n]];
        mask[ptr2] = true;
    }

    int nnz = 0;
    for (int i = 0; i<sza; i++){
        ptr1 = msh->mprts.parts[msh->mprts.parts_idx[part]+a[i]];
        for (j = msh->grph.row_idx[ptr1]; j<msh->grph.row_idx[ptr1+1]; j++){
            if (mask[msh->grph.col_idx[j]]){nnz++;}
        }
    }

    free(mask);
    return nnz;
}
