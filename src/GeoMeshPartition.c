#include "GeoMesh.h"

void stereo_Project_up(struct DoubleMatrix* points2d, struct DoubleMatrix* points3d);
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
    int nouter = (int) ceil(log((double) ntries-nlines+1) / lod(20.0));
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
        // find  C centerpoint points3d, csample

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