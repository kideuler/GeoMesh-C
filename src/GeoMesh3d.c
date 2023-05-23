#include "GeoMesh.h"

// global arrays
int Faces[4][3] = {{0,2,1},{0,1,3},{1,2,3},{2,0,3}};
int v2av[4][3] = {{1,3,2},{0,2,3},{3,1,0},{2,0,1}};
int v2f[4][3] = {{1,3,0},{0,2,1},{2,0,3},{3,1,2}};

// TODO:
// imrpove recursive tet find
// add finding ahf

struct Mesh GeoMesh_Delaunay_tet(struct DoubleMatrix *xs){
    int nv = xs->nrows;

    struct Mesh msh;
    msh.nelems = 0;
    int upper_bound = nv*8;
    msh.elems = IntMatrix_create(upper_bound,4);
    msh.sibhfs = IntMatrix_create(upper_bound,4);
    msh.on_boundary = (bool* ) malloc(upper_bound*sizeof(bool));
    memset(msh.on_boundary, false, upper_bound*sizeof(bool));
    msh.bwork = (bool* ) malloc(upper_bound*sizeof(bool));
    memset(msh.bwork, false, upper_bound*sizeof(bool));
    msh.stack = (int*) malloc(2*upper_bound*sizeof(int));

    // creating temporary buffers
    struct DoubleMatrix *coords = (struct DoubleMatrix *) malloc(sizeof(struct DoubleMatrix));
    coords->nrows = nv+4;
    coords->ncols = 3;
    coords->data = (double* ) malloc((nv+4)*3*sizeof(double));

    double ax = DoubleMatrix_min(xs,0); double ay = DoubleMatrix_min(xs,1); double az = DoubleMatrix_min(xs,2);
    double bx = DoubleMatrix_max(xs,0); double by = DoubleMatrix_max(xs,1); double bz = DoubleMatrix_max(xs,2);
    double d = Max(Max(bx-ax,by-ay),bz-az);

    for (int n = 0; n<nv; n++){
        coords->data[n*coords->ncols] = (xs->data[n*xs->ncols] - ax)/d + 1e-8*drand(-1.0,1.0);
        coords->data[n*coords->ncols+1] = (xs->data[n*xs->ncols+1] - ay)/d + 1e-8*drand(-1.0,1.0);
        coords->data[n*coords->ncols+2] = (xs->data[n*xs->ncols+2] - az)/d + 1e-8*drand(-1.0,1.0);
    }

    // sort points by binsort later

    // create big triangle
    coords->data[nv*3] = -100.0;
    coords->data[nv*3+1] = -100.0;
    coords->data[nv*3+2] = -100.0;

    coords->data[(nv+1)*3] = 100.0;
    coords->data[(nv+1)*3+1] = -100.0;
    coords->data[(nv+1)*3+2] = -100.0;

    coords->data[(nv+2)*3] = 0.0;
    coords->data[(nv+2)*3+1] = 100.0;
    coords->data[(nv+2)*3+2] = -100.0;

    coords->data[(nv+3)*3] = 0.0;
    coords->data[(nv+3)*3+1] = 0.0;
    coords->data[(nv+3)*3+2] = 100.0;

    msh.elems.data[0] = nv;
    msh.elems.data[1] = nv+1; 
    msh.elems.data[2] = nv+2;
    msh.elems.data[3] = nv+3;
    msh.sibhfs.data[0] = 0; 
    msh.sibhfs.data[1] = 0; 
    msh.sibhfs.data[2] = 0; 
    msh.sibhfs.data[3] = 0;

    // main loop inserting each node into the triangulation
    int enclosed_tet=-1;
    int vid;
    bool exitl, inside = true;
    msh.coords = *coords;
    msh.nelems = 1;
    msh.delete_elem = (bool*) malloc(upper_bound*sizeof(bool));
    memset(msh.delete_elem, false, upper_bound*sizeof(bool));
    for (int n = 0; n<nv; n++){
        vid = n;
        for (int i = msh.nelems-1; i>=0; i--){
            if (!msh.delete_elem[i]){
                enclosed_tet = i; break;
            }
        }
        // find enclosing tetrahedron
        double ps[3] = {msh.coords.data[vid*msh.coords.ncols], msh.coords.data[vid*msh.coords.ncols+1], msh.coords.data[vid*msh.coords.ncols+2]};
        inside = Mesh_find_enclosing_tet(&msh, &enclosed_tet, ps);
        
        if (!inside){ printf("no enclosing tetrahdron found\n");}

        // insert node into mesh
        Mesh_Bowyer_Watson_insertion3D(&msh, vid, enclosed_tet);
        printf("%d\n",n+1);
    }
    //free(bins);
    //free(order);


    // deleting any elements with big triangle points in them
    for (int i = 0; i<msh.nelems; i++){
        for (int j = 0; j<4; j++){
            if (msh.elems.data[i*4+j] >= nv){
                msh.delete_elem[i] = true;
            }
        }
    }
    Mesh_deleteElems3D(&msh);


    msh.coords = *xs;

    // freeing temporary buffers
    free(coords->data); free(coords);
    return msh;
}

void Bowyer_watson_recursive_tet_find(struct Mesh* msh, const double ps[3], int tet, int* stack_size);
void Mesh_Bowyer_Watson_insertion3D(struct Mesh* msh, int vid, int tet_start){
    int stack_size = -1;
    double xs[4][3], ps[3];
    ps[0] = msh->coords.data[3*vid]; ps[1] = msh->coords.data[3*vid+1]; ps[2] = msh->coords.data[3*vid+2]; 

    // find all tets whose circumsphere contains point
    Bowyer_watson_recursive_tet_find(msh, ps, tet_start, &stack_size);
    /*
    int v1,v2,v3,v4;
    for (int i = 0; i<msh->nelems; i++){
        if (!msh->delete_elem[i]){
            v1 = msh->elems.data[i*(msh->elems.ncols)];
            v2 = msh->elems.data[i*(msh->elems.ncols)+1];
            v3 = msh->elems.data[i*(msh->elems.ncols)+2];
            v4 = msh->elems.data[i*(msh->elems.ncols)+3];
            xs[0][0] = msh->coords.data[v1*msh->coords.ncols];
            xs[0][1] = msh->coords.data[v1*msh->coords.ncols+1];
            xs[0][2] = msh->coords.data[v1*msh->coords.ncols+2];

            xs[1][0] = msh->coords.data[v2*msh->coords.ncols];
            xs[1][1] = msh->coords.data[v2*msh->coords.ncols+1];
            xs[1][2] = msh->coords.data[v2*msh->coords.ncols+2];

            xs[2][0] = msh->coords.data[v3*msh->coords.ncols];
            xs[2][1] = msh->coords.data[v3*msh->coords.ncols+1];
            xs[2][2] = msh->coords.data[v3*msh->coords.ncols+2];

            xs[3][0] = msh->coords.data[v4*msh->coords.ncols];
            xs[3][1] = msh->coords.data[v4*msh->coords.ncols+1];
            xs[3][2] = msh->coords.data[v4*msh->coords.ncols+2];
            if (inside_circumtet(xs, ps)){
                stack_size++;
                msh->stack[stack_size] = i;
                msh->delete_elem[i] = true;
            }
        }
    }
    */

    // find all faces which form the convex polygon
    int eid, hfid, oppeid, opplid;
    int nelems_new = 0;
    for (int ii = 0; ii<=stack_size; ii++){
        eid = msh->stack[ii];
        for (int jj = 0; jj<4; jj++){
            hfid = msh->sibhfs.data[4*eid + jj];
            oppeid = (hfid==0)?0:hfid2eid(hfid)-1;
            if (hfid == 0 || !msh->delete_elem[oppeid]){
                msh->stack[4*nelems_new+stack_size+1] = msh->elems.data[4*eid + Faces[jj][0]];
                msh->stack[4*nelems_new+stack_size+2] = msh->elems.data[4*eid + Faces[jj][2]];
                msh->stack[4*nelems_new+stack_size+3] = msh->elems.data[4*eid + Faces[jj][1]];
                msh->stack[4*nelems_new+stack_size+4] = hfid;
                nelems_new++;
            }
        }
    }

    // reorganizing data for most efficient space storage
    int eids[nelems_new];

    int sz = (nelems_new<stack_size+1)?nelems_new:(stack_size+1);
    for (int ii = 0; ii<sz; ii++){
        eids[ii] = msh->stack[ii];
    }

    for (int ii = sz; ii<nelems_new;ii++){
        eids[ii] = msh->nelems + ii-sz;
    }
    msh->nelems+=(nelems_new-sz);
    
    // adding new elements to mesh
    for (int ii = 0; ii<nelems_new; ii++){
        eid = eids[ii];
        msh->delete_elem[eid] = false;
        msh->elems.data[4*eid] = msh->stack[4*ii+stack_size+1];
        msh->elems.data[4*eid+1] = msh->stack[4*ii+stack_size+2];
        msh->elems.data[4*eid+2] = msh->stack[4*ii+stack_size+3];
        msh->elems.data[4*eid+3] = vid;
        hfid = msh->stack[4*ii+stack_size+4];
        msh->sibhfs.data[4*eid] = hfid;
        msh->sibhfs.data[4*eid+1] = 0;
        msh->sibhfs.data[4*eid+2] = 0;
        msh->sibhfs.data[4*eid+3] = 0;
        if (hfid != 0){
            oppeid = hfid2eid(hfid)-1;
            opplid = hfid2lid(hfid)-1;
            msh->sibhfs.data[4*oppeid + opplid] = elids2hfid(eid+1,1);
        }
    }

    // update the ahf array
    Mesh_compute_AHF3D(msh);

    /*
    for (int ii = 0; ii<nelems_new; ii++){
        eid = eids[ii];
        v1 = msh->elems.data[eid*(msh->elems.ncols)];
        v2 = msh->elems.data[eid*(msh->elems.ncols)+1];
        v3 = msh->elems.data[eid*(msh->elems.ncols)+2];
        v4 = msh->elems.data[eid*(msh->elems.ncols)+3];
        
        xs[0][0] = msh->coords.data[v1*msh->coords.ncols];
        xs[0][1] = msh->coords.data[v1*msh->coords.ncols+1];
        xs[0][2] = msh->coords.data[v1*msh->coords.ncols+2];

        xs[1][0] = msh->coords.data[v2*msh->coords.ncols];
        xs[1][1] = msh->coords.data[v2*msh->coords.ncols+1];
        xs[1][2] = msh->coords.data[v2*msh->coords.ncols+2];

        xs[2][0] = msh->coords.data[v3*msh->coords.ncols];
        xs[2][1] = msh->coords.data[v3*msh->coords.ncols+1];
        xs[2][2] = msh->coords.data[v3*msh->coords.ncols+2];

        xs[3][0] = msh->coords.data[v4*msh->coords.ncols];
        xs[3][1] = msh->coords.data[v4*msh->coords.ncols+1];
        xs[3][2] = msh->coords.data[v4*msh->coords.ncols+2];
        printf("%d %d %d %d | %d %d %d %d | %f\n",msh->elems.data[4*eid], msh->elems.data[4*eid+1], msh->elems.data[4*eid+2],msh->elems.data[4*eid+3], \
        hfid2eid(msh->sibhfs.data[4*eid])-1, hfid2eid(msh->sibhfs.data[4*eid+1])-1, hfid2eid(msh->sibhfs.data[4*eid+2])-1,hfid2eid(msh->sibhfs.data[4*eid+3])-1, tet_vol(xs));
    }
    */
}

void Bowyer_watson_recursive_tet_find(struct Mesh* msh, const double ps[3], int tet, int* stack_size){
    double xs[4][3];
    int hfid,oppeid,v1,v2,v3,v4;
    v1 = msh->elems.data[tet*(msh->elems.ncols)];
    v2 = msh->elems.data[tet*(msh->elems.ncols)+1];
    v3 = msh->elems.data[tet*(msh->elems.ncols)+2];
    v4 = msh->elems.data[tet*(msh->elems.ncols)+3];
    xs[0][0] = msh->coords.data[v1*msh->coords.ncols];
    xs[0][1] = msh->coords.data[v1*msh->coords.ncols+1];
    xs[0][2] = msh->coords.data[v1*msh->coords.ncols+2];

    xs[1][0] = msh->coords.data[v2*msh->coords.ncols];
    xs[1][1] = msh->coords.data[v2*msh->coords.ncols+1];
    xs[1][2] = msh->coords.data[v2*msh->coords.ncols+2];

    xs[2][0] = msh->coords.data[v3*msh->coords.ncols];
    xs[2][1] = msh->coords.data[v3*msh->coords.ncols+1];
    xs[2][2] = msh->coords.data[v3*msh->coords.ncols+2];

    xs[3][0] = msh->coords.data[v4*msh->coords.ncols];
    xs[3][1] = msh->coords.data[v4*msh->coords.ncols+1];
    xs[3][2] = msh->coords.data[v4*msh->coords.ncols+2];

    if (inside_circumtet(xs,ps)){
        *stack_size = *stack_size + 1;
        msh->stack[(*stack_size)] = tet;
        msh->delete_elem[tet] = true;
        hfid = msh->sibhfs.data[4*tet];
        oppeid = hfid2eid(hfid)-1;
        if (hfid!=0 && !msh->delete_elem[oppeid]){
            Bowyer_watson_recursive_tet_find(msh, ps, oppeid, stack_size);
        }

        hfid = msh->sibhfs.data[4*tet+1];
        oppeid = hfid2eid(hfid)-1;
        if (hfid!=0 && !msh->delete_elem[oppeid]){
            Bowyer_watson_recursive_tet_find(msh, ps, oppeid, stack_size);
        }

        hfid = msh->sibhfs.data[4*tet+2];
        oppeid = hfid2eid(hfid)-1;
        if (hfid!=0 && !msh->delete_elem[oppeid]){
            Bowyer_watson_recursive_tet_find(msh, ps, oppeid, stack_size);
        }

        hfid = msh->sibhfs.data[4*tet+3];
        oppeid = hfid2eid(hfid)-1;
        if (hfid!=0 && !msh->delete_elem[oppeid]){
            Bowyer_watson_recursive_tet_find(msh, ps, oppeid, stack_size);
        }
    }
}

void Mesh_compute_AHF3D(struct Mesh* msh){
    int next[3] = {1,2,0};
    int prev[3] = {2,0,1};
    int nv = msh->coords.nrows;
    int nelems = msh->nelems;
    int nf = 4;
    memset(msh->sibhfs.data, 0, msh->sibhfs.nrows*4*sizeof(int));

    int* is_index = (int*) malloc((nv+1)*sizeof(int));
    memset(is_index, 0, (nv+1)*sizeof(int));
    int v,c,vn,ii,jj,kk;
    for (ii = 0; ii<nelems; ii++){
        if (!msh->delete_elem[ii]){
        for (jj = 0; jj<4; jj++){
            v = -1;
            for (kk=0; kk<3; kk++){
                if (msh->elems.data[4*ii + Faces[jj][kk]] > v){
                    v = msh->elems.data[4*ii + Faces[jj][kk]];
                }
            }
            is_index[v+1]++;
        }
        }
    }
    is_index[0] = 0;

    for(ii=0; ii<nv; ii++){is_index[ii+1] = is_index[ii]+is_index[ii+1];}

    int ne = nelems*nf;

    int* v2hf_fid = (int*) malloc(ne*sizeof(int));
    int* v2hf_leid = (int*) malloc(ne*sizeof(int));
    int* v2oe_v1 = (int*) malloc(ne*sizeof(int));
    int* v2oe_v2 = (int*) malloc(ne*sizeof(int));
    int av[3];
    for (ii = 0; ii<nelems; ii++){
        if (!msh->delete_elem[ii]){
        for (jj=0;jj<4;jj++){
            v = msh->elems.data[4*ii+jj];
            for(kk=0;kk<3;kk++){av[kk]=msh->elems.data[4*ii+v2av[jj][kk]];}
            for(kk=0;kk<3;kk++){
                if(v>av[kk] && v>av[next[kk]]){
                    v2oe_v1[is_index[v]] = av[kk];
                    v2oe_v2[is_index[v]] = av[next[kk]];
                    v2hf_fid[is_index[v]] = ii;
                    v2hf_leid[is_index[v]] = v2f[jj][kk];
                    is_index[v]++;
                }
            }
        }
        }
    }

    for(int ii=nv-1; ii>=0; ii--){ is_index[ii+1] = is_index[ii]; }
    is_index[0] = 0;
    int vs[3];
    int first_eid, first_lid, prev_eid, prev_lid, m_idx, nhfs, v1,v2;
    bool oriented = true;
    for (ii=0;ii<nelems;ii++){
        if (!msh->delete_elem[ii]){
        for (jj=0;jj<4;jj++){
            if (msh->sibhfs.data[4*ii+jj]>0){continue;}
            for(kk=0;kk<3;kk++){vs[kk] = msh->elems.data[4*ii + Faces[jj][kk]];}
            v = -1;
            for(kk=0;kk<3;kk++){
                if(vs[kk]>v){
                    v = vs[kk];
                    m_idx = kk;
                }
            }

            first_eid = ii; first_lid = jj;
            prev_eid = ii; prev_lid = jj;
            nhfs = 0;
            v1 = vs[prev[m_idx]]; v2 = vs[next[m_idx]];

            // search for sibling half-face in opposite direction
            for(int index = is_index[v]; index < is_index[v+1]; index++){
                if (v2oe_v1[index] == v1 && v2oe_v2[index] == v2){
                    msh->sibhfs.data[4*prev_eid + prev_lid] = elids2hfid(v2hf_fid[index]+1, v2hf_leid[index]+1);
                    prev_eid = v2hf_fid[index];
                    prev_lid = v2hf_leid[index];
                    nhfs++;
                }
            }

            // search for half-face in same orientation
            for(int index = is_index[v]; index < is_index[v+1]; index++){
                if (v2oe_v1[index] == v2 && v2oe_v2[index] == v1 && v2hf_fid[index] != ii){
                    msh->sibhfs.data[4*prev_eid + prev_lid] = elids2hfid(v2hf_fid[index]+1, v2hf_leid[index]+1);
                    prev_eid = v2hf_fid[index];
                    prev_lid = v2hf_leid[index];
                    nhfs++;
                    oriented = false;
                }
            }

            // close up cycle
            if (prev_eid != first_eid){
                msh->sibhfs.data[4*prev_eid + prev_lid] = elids2hfid(first_eid+1,first_lid+1);
                nhfs++;
            }
        }
        }
    }

    if (!oriented){
        printf("mesh is not oriented\n");
    }


    free(is_index);
    free(v2hf_fid);
    free(v2hf_leid);
    free(v2oe_v1);
    free(v2oe_v2);
}

void Mesh_deleteElems3D(struct Mesh* msh){
    
    int nelems=0;
    for (int i = 0; i<msh->nelems; i++){
        if (!msh->delete_elem[i]){nelems++;}
    }
    int* elems_data  = (int*) malloc(4*nelems*sizeof(int));
    int* sibhfs_data  = (int*) malloc(4*nelems*sizeof(int));
    int* idx = (int*) malloc(msh->nelems*sizeof(int));
    int* idx_rev = (int*) malloc(msh->nelems*sizeof(int));
    memset(idx_rev,-1,msh->nelems*sizeof(int));

    int k = 0;
    for (int i = 0; i<msh->nelems; i++){
        if (!msh->delete_elem[i]) {
            for (int j = 0; j<4; j++){
                elems_data[4*k + j] = msh->elems.data[4*i + j];
                sibhfs_data[4*k+j] = msh->sibhfs.data[4*i+j];
            }
            msh->delete_elem[k] = false;
            idx[k] = i;
            k++;
        } else {
            msh->delete_elem[i] = false;
        }
    }
    msh->nelems = nelems;

    for (int i = 0; i<msh->nelems; i++){
        idx_rev[idx[i]] = i;
    }

    // fix sibhfs
    int nside;
    int hfid, eid, lid;
    for (int i = 0; i<nelems; i++){
        nside = 0;
        for (int j = 0; j < 4; j++){
            hfid = sibhfs_data[4*i+j];
            if (!hfid == 0){
                eid = hfid2eid(hfid);
                lid = hfid2lid(hfid);
                if (!msh->delete_elem[eid-1]){
                    if (idx_rev[eid-1] == -1){
                        sibhfs_data[4*i+j] = 0;
                    } else {
                        sibhfs_data[4*i+j] = elids2hfid(idx_rev[eid-1]+1,lid);
                    }
                } else {
                    sibhfs_data[4*i+j] = 0;
                }
                nside++;
            } else {
                sibhfs_data[4*i+j] = 0;
            }
        }
        msh->delete_elem[i] = false;
    }


    free(idx); free(idx_rev);
    free(msh->elems.data);
    msh->elems.data = elems_data;
    msh->elems.nrows = nelems;
    free(msh->sibhfs.data);
    msh->sibhfs.nrows = nelems;
    msh->sibhfs.data = sibhfs_data;
}

bool Mesh_find_enclosing_tet_noAHF(struct Mesh* msh, int* tet, double ps[3]){
    int v1,v2,v3,v4,i;
    double xs[4][3];
    for (i = 0; i<msh->nelems; i++){
        if (!msh->delete_elem[i]){
            v1 = msh->elems.data[i*(msh->elems.ncols)];
            v2 = msh->elems.data[i*(msh->elems.ncols)+1];
            v3 = msh->elems.data[i*(msh->elems.ncols)+2];
            v4 = msh->elems.data[i*(msh->elems.ncols)+3];
            xs[0][0] = msh->coords.data[v1*msh->coords.ncols];
            xs[0][1] = msh->coords.data[v1*msh->coords.ncols+1];
            xs[0][2] = msh->coords.data[v1*msh->coords.ncols+2];

            xs[1][0] = msh->coords.data[v2*msh->coords.ncols];
            xs[1][1] = msh->coords.data[v2*msh->coords.ncols+1];
            xs[1][2] = msh->coords.data[v2*msh->coords.ncols+2];

            xs[2][0] = msh->coords.data[v3*msh->coords.ncols];
            xs[2][1] = msh->coords.data[v3*msh->coords.ncols+1];
            xs[2][2] = msh->coords.data[v3*msh->coords.ncols+2];

            xs[3][0] = msh->coords.data[v4*msh->coords.ncols];
            xs[3][1] = msh->coords.data[v4*msh->coords.ncols+1];
            xs[3][2] = msh->coords.data[v4*msh->coords.ncols+2];

            if(inside_tet(xs,ps)==0){
                *tet = i;
                return true;
            }
        }
    }
    return false;
}

bool Mesh_find_enclosing_tet(struct Mesh* msh, int* tet, double ps[3]){
    int v1,v2,v3,v4,i,hfid;
    bool stop;
    int iters = 0;
    i = 0;
    stop = false;
    double xs[4][3] = {{0,0,0},{0,0,0},{0,0,0},{0,0,0}};
    int IT;
    while (!stop){
        v1 = msh->elems.data[4*(*tet)];
        v2 = msh->elems.data[4*(*tet)+1];
        v3 = msh->elems.data[4*(*tet)+2];
        v4 = msh->elems.data[4*(*tet)+3];
        if (msh->delete_elem[*tet]){
            printf("deleted tet passed ran through in find_enclosing_tet, should not be\n");

        }
        xs[0][0] = msh->coords.data[v1*msh->coords.ncols];
        xs[0][1] = msh->coords.data[v1*msh->coords.ncols+1];
        xs[0][2] = msh->coords.data[v1*msh->coords.ncols+2];

        xs[1][0] = msh->coords.data[v2*msh->coords.ncols];
        xs[1][1] = msh->coords.data[v2*msh->coords.ncols+1];
        xs[1][2] = msh->coords.data[v2*msh->coords.ncols+2];

        xs[2][0] = msh->coords.data[v3*msh->coords.ncols];
        xs[2][1] = msh->coords.data[v3*msh->coords.ncols+1];
        xs[2][2] = msh->coords.data[v3*msh->coords.ncols+2];

        xs[3][0] = msh->coords.data[v4*msh->coords.ncols];
        xs[3][1] = msh->coords.data[v4*msh->coords.ncols+1];
        xs[3][2] = msh->coords.data[v4*msh->coords.ncols+2];

        IT = inside_tet(xs,ps);
        if (IT==0){
            stop = true;
            return true;
        } else if (IT==1){
            hfid = msh->sibhfs.data[4*(*tet)];
            if (hfid != 0){
                *tet = hfid2eid(hfid)-1;
            } else {
                stop = true;
                *tet = elids2hfid(*tet+1,1);
                return false;
            }
        } else if (IT==2){
            hfid = msh->sibhfs.data[4*(*tet)+1];
            if (hfid != 0){
                *tet = hfid2eid(hfid)-1;
            } else {
                stop = true;
                *tet = elids2hfid(*tet+1,1);
                return false;
            }
        } else if (IT==3){
           hfid = msh->sibhfs.data[4*(*tet)+2];
            if (hfid != 0){
                *tet = hfid2eid(hfid)-1;
            } else {
                stop = true;
                *tet = elids2hfid(*tet+1,1);
                return false;
            }
        } else if (IT==4){
            hfid = msh->sibhfs.data[4*(*tet)+3];
            if (hfid != 0){
                *tet = hfid2eid(hfid)-1;
            } else {
                stop = true;
                *tet = elids2hfid(*tet+1,1);
                return false;
            }
        }
        iters++;
    }
    *tet = -1;
    return false;
    
}

void normal3d(double* N, const double xs[4][3], int facet){
    double a[3], b[3];
    a[0] = xs[Faces[facet][1]][0] - xs[Faces[facet][0]][0];
    a[1] = xs[Faces[facet][1]][1] - xs[Faces[facet][0]][1];
    a[2] = xs[Faces[facet][1]][2] - xs[Faces[facet][0]][2];
    b[0] = xs[Faces[facet][2]][0] - xs[Faces[facet][0]][0];
    b[1] = xs[Faces[facet][2]][1] - xs[Faces[facet][0]][1];
    b[2] = xs[Faces[facet][2]][2] - xs[Faces[facet][0]][2];

    N[0] = a[1]*b[2] - a[2]*b[1];
    N[1] = a[2]*b[0] - a[0]*b[2];
    N[2] = a[0]*b[1] - a[1]*b[0];
}

int inside_tet(const double xs[4][3], const double ps[3]){

    double* a = (double*) malloc(3*sizeof(double));
    double* p = (double*) malloc(3*sizeof(double));
    double* N = (double*) malloc(3*sizeof(double));

    normal3d(N,xs,0);
    a[0] = xs[3][0]-xs[Faces[0][0]][0];
    a[1] = xs[3][1]-xs[Faces[0][0]][1];
    a[2] = xs[3][2]-xs[Faces[0][0]][2];
    double D0 = N[0]*a[0]+N[1]*a[1]+N[2]*a[2];
    p[0] = ps[0]-xs[Faces[0][0]][0];
    p[1] = ps[1]-xs[Faces[0][0]][1];
    p[2] = ps[2]-xs[Faces[0][0]][2];

    double D1 = N[0]*p[0]+N[1]*p[1]+N[2]*p[2];
    bool S1 = D1 < 0;

    normal3d(N,xs,1);
    a[0] = xs[2][0]-xs[Faces[1][0]][0];
    a[1] = xs[2][1]-xs[Faces[1][0]][1];
    a[2] = xs[2][2]-xs[Faces[1][0]][2];
    D0 = N[0]*a[0]+N[1]*a[1]+N[2]*a[2];
    p[0] = ps[0]-xs[Faces[1][0]][0];
    p[1] = ps[1]-xs[Faces[1][0]][1];
    p[2] = ps[2]-xs[Faces[1][0]][2];
    double D2 = N[0]*p[0]+N[1]*p[1]+N[2]*p[2];
    bool S2 = D2 < 0;

    normal3d(N,xs,2);
    a[0] = xs[0][0]-xs[Faces[2][0]][0];
    a[1] = xs[0][1]-xs[Faces[2][0]][1];
    a[2] = xs[0][2]-xs[Faces[2][0]][2];
    D0 = N[0]*a[0]+N[1]*a[1]+N[2]*a[2];
    p[0] = ps[0]-xs[Faces[2][0]][0];
    p[1] = ps[1]-xs[Faces[2][0]][1];
    p[2] = ps[2]-xs[Faces[2][0]][2];
    double D3 = N[0]*p[0]+N[1]*p[1]+N[2]*p[2];
    bool S3 = D3 < 0;

    normal3d(N,xs,3);
    a[0] = xs[1][0]-xs[Faces[3][0]][0];
    a[1] = xs[1][1]-xs[Faces[3][0]][1];
    a[2] = xs[1][2]-xs[Faces[3][0]][2];
    D0 = N[0]*a[0]+N[1]*a[1]+N[2]*a[2];
    p[0] = ps[0]-xs[Faces[3][0]][0];
    p[1] = ps[1]-xs[Faces[3][0]][1];
    p[2] = ps[2]-xs[Faces[3][0]][2];
    double D4 = N[0]*p[0]+N[1]*p[1]+N[2]*p[2];
    bool S4 = D4 < 0;
    free(a); free(p); free(N);
    if (S1 && S2 && S3 && S4){
        return 0;
    } else if(D1>=0 && (D1>=D2) && (D1>=D3) && (D1>=D4)){
        return 1;
    } else if(D2>=0 && (D2>=D1) && (D2>=D3) && (D2>=D4)){
        return 2;
    } else if(D3>=0 && (D3>=D1) && (D3>=D2) && (D3>=D4)){
        return 3;
    } else if(D4>=0 && (D4>=D1) && (D4>=D2) && (D4>=D3)){
        return 4;
    } else {
        printf("error in inside_tet returning null");
        assert(false);
        return -1;
    }
}

bool inside_circumtet(const double xs[4][3], const double ps[3]){
    double* C = (double*) malloc(3*sizeof(double));
    circumcenter_tet(xs,C);
    double R = (xs[0][0]-C[0])*(xs[0][0]-C[0]) + (xs[0][1]-C[1])*(xs[0][1]-C[1]) + (xs[0][2]-C[2])*(xs[0][2]-C[2]);
    bool D = ((ps[0]-C[0])*(ps[0]-C[0]) + (ps[1]-C[1])*(ps[1]-C[1]) + (ps[2]-C[2])*(ps[2]-C[2])) <= R;
    free(C);
    return D;
}

void circumcenter_tet(const double xs[4][3], double* C){
    double denominator;
 
    // Use coordinates relative to point `a' of the tetrahedron.
 
    // ba = b - a
    double ba_x = xs[1][0] - xs[0][0];
    double ba_y = xs[1][1] - xs[0][1];
    double ba_z = xs[1][2] - xs[0][2];
 
    // ca = c - a
    double ca_x = xs[2][0] - xs[0][0];
    double ca_y = xs[2][1] - xs[0][1];
    double ca_z = xs[2][2] - xs[0][2];
 
    // da = d - a
    double da_x = xs[3][0] - xs[0][0];
    double da_y = xs[3][1] - xs[0][1];
    double da_z = xs[3][2] - xs[0][2];
 
    // Squares of lengths of the edges incident to `a'.
    double len_ba = ba_x * ba_x + ba_y * ba_y + ba_z * ba_z;
    double len_ca = ca_x * ca_x + ca_y * ca_y + ca_z * ca_z;
    double len_da = da_x * da_x + da_y * da_y + da_z * da_z;
 
    // Cross products of these edges.
 
    // c cross d
    double cross_cd_x = ca_y * da_z - da_y * ca_z;
    double cross_cd_y = ca_z * da_x - da_z * ca_x;
    double cross_cd_z = ca_x * da_y - da_x * ca_y;
 
    // d cross b
    double cross_db_x = da_y * ba_z - ba_y * da_z;
    double cross_db_y = da_z * ba_x - ba_z * da_x;
    double cross_db_z = da_x * ba_y - ba_x * da_y;
 
    // b cross c
    double cross_bc_x = ba_y * ca_z - ca_y * ba_z;
    double cross_bc_y = ba_z * ca_x - ca_z * ba_x;
    double cross_bc_z = ba_x * ca_y - ca_x * ba_y;
 
    // Calculate the denominator of the formula.
    denominator = 0.5 / (ba_x * cross_cd_x + ba_y * cross_cd_y + ba_z * cross_cd_z);
 
    // Calculate offset (from `a') of circumcenter.
    double circ_x = (len_ba * cross_cd_x + len_ca * cross_db_x + len_da * cross_bc_x) * denominator;
    double circ_y = (len_ba * cross_cd_y + len_ca * cross_db_y + len_da * cross_bc_y) * denominator;
    double circ_z = (len_ba * cross_cd_z + len_ca * cross_db_z + len_da * cross_bc_z) * denominator;
    
    C[0] = circ_x + xs[0][0];
    C[1] = circ_y + xs[0][1];
    C[2] = circ_z + xs[0][2];
    return;
}

static double dphi[4][3] = {{-1.0,-1.0,-1.0},{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
double tet_vol(const double xs[4][3]){
    double J[3][3];

    for (int ii = 0; ii<3; ii++){
        for (int jj=0; jj<3;jj++){
            J[ii][jj] = 0.0;
            for (int kk=0; kk<4; kk++){
                J[ii][jj] += xs[kk][ii]*dphi[kk][jj];
            }
        }
    }

    double detJ = J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1]) - \
    J[0][1]*(J[1][0]*J[2][2] - J[1][2]*J[2][0]) + \
    J[0][2]*(J[1][0]*J[2][1] - J[2][0]*J[1][1]);

    return detJ;
}

bool check_jacobians_tet(struct Mesh* msh){
    bool check = true;
    int v1,v2,v3,v4,eid;
    double xs[4][3], vol;
    for (int ii = 0; ii<msh->nelems; ii++){
        eid = ii;
        v1 = msh->elems.data[eid*(msh->elems.ncols)];
        v2 = msh->elems.data[eid*(msh->elems.ncols)+1];
        v3 = msh->elems.data[eid*(msh->elems.ncols)+2];
        v4 = msh->elems.data[eid*(msh->elems.ncols)+3];
        
        xs[0][0] = msh->coords.data[v1*msh->coords.ncols];
        xs[0][1] = msh->coords.data[v1*msh->coords.ncols+1];
        xs[0][2] = msh->coords.data[v1*msh->coords.ncols+2];

        xs[1][0] = msh->coords.data[v2*msh->coords.ncols];
        xs[1][1] = msh->coords.data[v2*msh->coords.ncols+1];
        xs[1][2] = msh->coords.data[v2*msh->coords.ncols+2];

        xs[2][0] = msh->coords.data[v3*msh->coords.ncols];
        xs[2][1] = msh->coords.data[v3*msh->coords.ncols+1];
        xs[2][2] = msh->coords.data[v3*msh->coords.ncols+2];

        xs[3][0] = msh->coords.data[v4*msh->coords.ncols];
        xs[3][1] = msh->coords.data[v4*msh->coords.ncols+1];
        xs[3][2] = msh->coords.data[v4*msh->coords.ncols+2];

        vol = tet_vol(xs);
        if (vol < 0){
            printf("negative jacobian at eid: %d volume: %f\n",ii,vol);
            check = false;
        }
    }
}

void Mesh2vtk_tet(struct Mesh* msh){
    FILE *fid;
    fid = fopen("test.vtk","w");
    fprintf(fid,"# vtk DataFile Version 3.0\n");
    fprintf(fid,"This file was written using writevtk_unstr.m\n");
    fprintf(fid,"ASCII\n");

    int ndims = msh->coords.ncols;
    int nv = msh->coords.nrows;
    // header for points
    fprintf(fid, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(fid, "POINTS %i double",nv);

    // write out vertices
    if (ndims == 2){
        for (int i=0; i<nv;i++){
            fprintf(fid,"\n%g %g %g",msh->coords.data[2*i],msh->coords.data[2*i+1],0.0);
        }
    } else {
        for (int i=0; i<nv;i++){
            fprintf(fid,"\n%g %g %g",msh->coords.data[3*i],msh->coords.data[3*i+1],msh->coords.data[3*i+2]);
        }
    }

    // write out connectivity header
    int nelems = msh->nelems;
    fprintf(fid,"\n\nCELLS %i %i", nelems, 5*nelems);
    for (int i = 0; i<nelems; i++){
        fprintf(fid,"\n%d %d %d %d %d",4,msh->elems.data[4*i],msh->elems.data[4*i+1],msh->elems.data[4*i+2], msh->elems.data[4*i+3]);
    }

    // write out cell types
    fprintf(fid, "\n\nCELL_TYPES %i", nelems);
    for (int i = 0; i<nelems; i++){
        fprintf(fid,"\n%i",10);
    }

    fclose(fid);
}