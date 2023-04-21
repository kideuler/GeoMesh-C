#include "GeoMesh.h"

int elids2hfid(int eid, int lid){
    return ((eid << 8) + lid - 1);
}

int hfid2eid(int hfid){
    return (hfid >> 8);
}

int hfid2lid(int hfid){
    return (hfid&255)+1;
}

void Mesh_print(struct Mesh *msh){
    printf("---ELEMENTS---\n");
    IntMatrix_print(&msh->elems);
    printf("---COORDINATES---\n");
    DoubleMatrix_print(&msh->coords);
    printf("---SIBHFS---\n");
    IntMatrix_print(&msh->sibhfs);
}

// Constrained Delaunay triangulation driver function
struct Mesh GeoMesh_ConstrainedDelaunay(struct IntMatrix *segments, struct DoubleMatrix *xs){
    struct Mesh msh = GeoMesh_Delaunay(xs);
    int maxne=10;
    Mesh_compute_OneringElements(&msh,maxne);
    maxne = msh.stncl.maxne;

    int* facets = (int*) malloc(2*msh.nelems*sizeof(int));
    msh.bwork = (bool*) malloc(msh.nelems*sizeof(bool));
    msh.delete_elem = (bool*) malloc(msh.nelems*sizeof(bool));

    int v1,v2,hvid,kk,eid,nid,nf,hfid,lid,oppeid,opplid;
    bool inray,connected,convex;
    nf = 0;
    for (int seg = 0; seg<segments->nrows; seg++){
        v1 = segments->data[2*seg]; v2 = segments->data[2*seg+1];
        connected = false;

        while (!connected){ // while segment not done
            inray = false;
            convex = false;
            hvid = 1;
            kk = 0;
            // finds segment or line that obstructs segments
            while (!inray && !connected && kk < maxne){
                hvid = msh.stncl.hvids.data[v1*maxne + kk];
                if (hvid <= 0){
                    break;
                }

                // check if either edge of the triangle connects to the targeted vertex
                eid = hfid2eid(hvid)-1;
                nid = hfid2lid(hvid)-1;
                if (msh.elems.data[3*eid + (nid+1)%3] == v2){
                    connected = true;
                    hfid = msh.sibhfs.data[3*eid+nid];
                    if (hfid > 0){
                        facets[2*nf] = eid;
                        facets[2*nf+1] = nid;
                        msh.bwork[eid] = true;
                        nf++;
                    }
                    break;
                }
                if (msh.elems.data[3*eid + (nid+2)%3] == v2){
                    connected = true;
                    hfid = msh.sibhfs.data[3*eid+(nid+2)%3];
                    if (hfid > 0){
                        facets[2*nf] = hfid2eid(hfid)-1;
                        facets[2*nf+1] = hfid2lid(hfid)-1;
                        msh.bwork[ hfid2eid(hfid)-1] = true;
                        nf++;
                    }
                    break;
                }

                // check if the opposite edge is obstructing the desired segment
                if (Line_cross(msh.coords.data[2*msh.elems.data[3*eid + (nid+1)%3]],msh.coords.data[2*msh.elems.data[3*eid + (nid+1)%3]+1], \
                msh.coords.data[2*msh.elems.data[3*eid + (nid+2)%3]],msh.coords.data[2*msh.elems.data[3*eid + (nid+2)%3]+1], \
                msh.coords.data[2*v1],msh.coords.data[2*v1+1], msh.coords.data[2*v2],msh.coords.data[2*v2+1])){
                    inray = true;
                    break;
                }
                kk++;
            }

            if (inray) { // flip segment and correct onering if convex quad
                // IMPLEMENT
            }
        }
    }

    // delete elements on opposite side of boundary
    if (true){
        for (int i=0; i<nf; i++){
            eid = facets[2*i];
            lid = facets[2*i+1];
            hfid = msh.sibhfs.data[3*eid+lid];
            Recursive_findDelete(&msh, hfid);
        }
    }
    Mesh_deleteElems(&msh);

    free(facets);

    return msh;
}

// Delaunay triangulation driver function
struct Mesh GeoMesh_Delaunay(struct DoubleMatrix *xs){
    int nv = xs->nrows;

    struct Mesh msh;
    msh.nelems = 0;
    int upper_bound = nv*4;
    msh.elems = IntMatrix_create(upper_bound,3);
    msh.sibhfs = IntMatrix_create(upper_bound,3);

    // creating temporary buffers
    struct DoubleMatrix *coords = (struct DoubleMatrix *) malloc(sizeof(struct DoubleMatrix));
    coords->nrows = nv+3;
    coords->ncols = 2;
    coords->data = (double* ) malloc((nv+3)*2*sizeof(double));

    double ax = DoubleMatrix_min(xs,0); double ay = DoubleMatrix_min(xs,1);
    double bx = DoubleMatrix_max(xs,0); double by = DoubleMatrix_max(xs,1);
    double d = bx-ax > by-bx ? bx-ax : by-ay;

    for (int n = 0; n<nv; n++){
        coords->data[n*coords->ncols] = (xs->data[n*xs->ncols] - ax)/d;
        coords->data[n*coords->ncols+1] = (xs->data[n*xs->ncols+1] - ay)/d;
    }

    // sort points by proximity using binsort
    int nbin = (int)ceil(pow(nv,0.5));
    int* bins = (int*)malloc(nv*sizeof(int));
    int* order = (int*)malloc(nv*sizeof(int));
    for (int n = 0; n<nv; n++){
        int p = (int)(coords->data[2*n]*nbin*0.999);
        int q = (int)(coords->data[2*n+1]*nbin*0.999);
        if (p%2){
            bins[n] = (p+1)*nbin-q;
        } else {
            bins[n] = p*nbin+q+1;
        } 
        bins[n] = (p%2) ? (p+1)*nbin-q : p*nbin+q+1;
        order[n] = n;
    }

    int key;
    int temp;
    int i,j;
    for (i = 1; i<nv; i++){
        key = bins[i];
        temp = order[i];
        j = i-1;
        while(j>=0 && bins[j]>key){
            bins[j+1] = bins[j];
            order[j+1] = order[j];
            j--;
        }
        bins[j+1] = key;
        order[j+1] = temp;
    }

    // create big triangle
    coords->data[nv*coords->ncols] = -100.0;
    coords->data[nv*coords->ncols+1] = -100.0;
    coords->data[(nv+1)*coords->ncols] = 100.0;
    coords->data[(nv+1)*coords->ncols+1] = -100.0;
    coords->data[(nv+2)*coords->ncols] = 0.0;
    coords->data[(nv+2)*coords->ncols+1] = 100.0;

    msh.elems.data[0] = nv;
    msh.elems.data[1] = nv+1; 
    msh.elems.data[2] = nv+2;
    msh.sibhfs.data[0] = 0; msh.sibhfs.data[1] = 0; msh.sibhfs.data[2] = 0;

    // main loop inserting each node into the triangulation
    int enclosed_tri=-1;
    int vid;
    bool exitl, inside;
    msh.coords = *coords;
    msh.nelems = 1;
    msh.delete_elem = (bool *) malloc(upper_bound*sizeof(bool));
    for (int n = 0; n<nv; n++){
        vid = order[n];
        enclosed_tri = msh.nelems-1;
        // find enclosing triangle
        double ps[2] = {msh.coords.data[vid*msh.coords.ncols],msh.coords.data[vid*msh.coords.ncols+1]};
        inside = Mesh_find_enclosing_tri(&msh, &enclosed_tri, ps);

        if (!inside){ printf("no enclosing triangle found\n");}

        // insert node into triangulation using sloans flip insertion algorithm
        Mesh_flip_insertion(&msh, &vid, enclosed_tri);
    }
    free(bins);
    free(order);


    // deleting any elements with big triangle points in them
    for (int i = 0; i<msh.nelems; i++){
        for (int j = 0; j<3; j++){
            if (msh.elems.data[i*3+j] >= nv){
                msh.delete_elem[i] = true;
            }
        }
    }


    Mesh_deleteElems(&msh);


    msh.coords = *xs;

    // freeing temporary buffers
    free(coords->data); free(coords);
    return msh;
}

void Mesh_compute_OneringElements(struct Mesh* msh, int maxne){
    int nv = msh->coords.nrows;
    msh->stncl.maxne = maxne;
    msh->stncl.hvids = IntMatrix_create(nv,maxne);
    memset(msh->stncl.hvids.data, -1, nv*maxne*sizeof(int));
    msh->stncl.nhvids = (int*)malloc(nv*sizeof(int));
    memset(msh->stncl.nhvids, 0, nv*sizeof(int));

    int v;
    for (int i = 0; i<msh->nelems; i++){
        for (int j = 0; j<3; j++){
            v = msh->elems.data[3*i + j];
            msh->stncl.nhvids[v]++;
            if ( msh->stncl.nhvids[v] >= maxne){
                free(msh->stncl.nhvids);
                free(msh->stncl.hvids.data);
                printf("Mesh_compute_OneringElements: buffers too small, enlarging buffers and rerunning\n");
                Mesh_compute_OneringElements(msh, 2*maxne);
                return;
            }
            msh->stncl.hvids.data[maxne*v + msh->stncl.nhvids[v]-1] = elids2hfid(i+1,j+1);
        }
    }

    return;
}

void Mesh_deleteElems(struct Mesh* msh){

    int nelems=0;
    for (int i = 0; i<msh->nelems; i++){
        if (!msh->delete_elem[i]){nelems++;}
    }

    int* elems_data  = (int*)malloc(3*nelems*sizeof(int));
    int* sibhfs_data  = (int*)malloc(3*nelems*sizeof(int));
    int* idx = (int*)malloc(msh->nelems*sizeof(int));
    int* idx_rev = (int*)malloc(msh->nelems*sizeof(int));
    memset(idx_rev,-1,nelems*sizeof(int));

    int k = 0;
    for (int i = 0; i<msh->nelems; i++){
        if (!msh->delete_elem[i]) {
            for (int j = 0; j<3; j++){
                elems_data[3*k] = msh->elems.data[3*i];
                elems_data[3*k+1] = msh->elems.data[3*i+1];
                elems_data[3*k+2] = msh->elems.data[3*i+2];
                sibhfs_data[3*k] = msh->sibhfs.data[3*i];
                sibhfs_data[3*k+1] = msh->sibhfs.data[3*i+1];
                sibhfs_data[3*k+2] = msh->sibhfs.data[3*i+2];
            }
            msh->delete_elem[k] = false;
            idx[k] = i;
            k++;
        } else {
            msh->delete_elem[i] = false;
        }
    }

    for (int i = 0; i<nelems; i++){
        idx_rev[idx[i]] = i;
    }

    // fix sibhfs
    int nside;
    int hfid, eid, lid;
    for (int i = 0; i<nelems; i++){
        nside = 0;
        for (int j = 0; j < 3; j++){
            hfid = sibhfs_data[3*i+j];
            if (!hfid == 0){
                eid = hfid2eid(hfid);
                lid = hfid2lid(hfid);
                if (!msh->delete_elem[eid-1]){
                    if (idx_rev[eid-1] == -1){
                        sibhfs_data[3*i+j] = 0;
                    } else {
                        sibhfs_data[3*i+j] = elids2hfid(idx_rev[eid-1]+1,lid);
                    }
                } else {
                    sibhfs_data[3*i+j] = 0;
                }
                nside++;
            } else {
                sibhfs_data[3*i+j] = 0;
            }
        }
        msh->delete_elem[i] = false;
    }


    free(idx); free(idx_rev);
    free(msh->elems.data);
    msh->elems.data = elems_data;
    free(msh->sibhfs.data);
    msh->sibhfs.data = sibhfs_data;
    msh->nelems = nelems;
}

void Recursive_findDelete(struct Mesh* msh, int hfid){
    int eid = hfid2eid(hfid)-1;
    int lid = hfid2lid(hfid)-1;
    if (!msh->delete_elem[eid]){
        msh->delete_elem[eid] = true;
        int hfid1 = msh->sibhfs.data[3*eid+(lid+1)%3];
        if (hfid1 != 0){
            if (!msh->bwork[hfid2eid(hfid1)-1]){
                Recursive_findDelete(msh, hfid1);
            }
        }
        int hfid2 = msh->sibhfs.data[3*eid+(lid+2)%3];
        if (hfid2 != 0){
            if (!msh->bwork[hfid2eid(hfid2)-1]){
                Recursive_findDelete(msh, hfid2);
            }
        }
    } 
    return;
}

void Mesh_flip_insertion(struct Mesh* msh, int* vid, int tri_start){
    msh->delete_elem[tri_start] = true;
    int hfid,eid,lid;

    int tri[3] = {msh->elems.data[tri_start*msh->elems.ncols],msh->elems.data[tri_start*msh->elems.ncols+1],msh->elems.data[tri_start*msh->elems.ncols+2]};
    int sib[3] = {msh->sibhfs.data[tri_start*msh->sibhfs.ncols],msh->sibhfs.data[tri_start*msh->sibhfs.ncols+1],msh->sibhfs.data[tri_start*msh->sibhfs.ncols+2]};

    int* stack = (int*) malloc(msh->elems.nrows*sizeof(int));
    int stack_size = -1;
    int eids[3] = {msh->nelems, msh->nelems+1, msh->nelems+2};
    // splitting triangles and adding them to the stack
    for (int i = 0; i<3; i++){
        msh->elems.data[eids[i]*msh->elems.ncols] = *vid;
        msh->elems.data[eids[i]*msh->elems.ncols+1] = tri[i];
        msh->elems.data[eids[i]*msh->elems.ncols+2] = tri[(i+1)%3];
        hfid = sib[i];
        msh->sibhfs.data[eids[i]*msh->sibhfs.ncols] = elids2hfid(eids[(i+2)%3]+1 ,3);
        msh->sibhfs.data[eids[i]*msh->sibhfs.ncols+1] = hfid;
        msh->sibhfs.data[eids[i]*msh->sibhfs.ncols+2] = elids2hfid(eids[(i+1)%3]+1,1);
        if (hfid2eid(hfid) > 0){
            msh->sibhfs.data[(hfid2eid(hfid)-1)*msh->sibhfs.ncols + hfid2lid(hfid)-1] = elids2hfid(eids[i]+1, 2);
            stack[stack_size+1] = elids2hfid(eids[i]+1, 2);
            stack_size++;
        }
    }

    msh->nelems+=3;
    double xs[3][2], ps[2];
    int oppeid, opplid;
    while (stack_size > -1){
        hfid = stack[stack_size];
        stack_size--;
        eid = hfid2eid(hfid)-1;
        lid = hfid2lid(hfid)-1;
        oppeid = hfid2eid(msh->sibhfs.data[eid*3 + lid])-1;
        opplid = hfid2lid(msh->sibhfs.data[eid*3 + lid])-1;
        xs[0][0] = msh->coords.data[2*msh->elems.data[3*oppeid]];
        xs[0][1] = msh->coords.data[2*msh->elems.data[3*oppeid]+1];
        xs[1][0] = msh->coords.data[2*msh->elems.data[3*oppeid+1]];
        xs[1][1] = msh->coords.data[2*msh->elems.data[3*oppeid+1]+1];
        xs[2][0] = msh->coords.data[2*msh->elems.data[3*oppeid+2]];
        xs[2][1] = msh->coords.data[2*msh->elems.data[3*oppeid+2]+1];
        ps[0] = msh->coords.data[2*msh->elems.data[3*eid]];
        ps[1] = msh->coords.data[2*msh->elems.data[3*eid]+1];
        if (inside_circumtri(xs,ps)){
            // flip edge
            flip_edge(msh,eid,1);
            
            hfid = msh->sibhfs.data[3*oppeid+1];
            if (hfid2eid(hfid) > 0){
                stack[stack_size+1] = elids2hfid(oppeid+1,2);
                stack_size++;
            }
            hfid = msh->sibhfs.data[3*eid+1];
            if (hfid2eid(hfid) > 0){
                stack[stack_size+1] = elids2hfid(eid+1,2);
                stack_size++;
            }
            
        }
    }


    free(stack);
    return;
}

bool Mesh_find_enclosing_tri(struct Mesh* msh, int* tri, double ps[2]){
    int v1,v2,v3,i,hfid;
    bool stop;
    int iters = 0;
    i = 0;
    stop = false;
    double xs[3][2] = {{0,0},{0,0},{0,0}};
    while (!stop && iters<10000){
        v1 = msh->elems.data[(*tri)*(msh->elems.ncols)];
        v2 = msh->elems.data[(*tri)*(msh->elems.ncols)+1];
        v3 = msh->elems.data[(*tri)*(msh->elems.ncols)+2];
        if (msh->delete_elem[*tri]){ printf("delete tri passed ran through Mesh_find_enclosing_tri\n");}

        xs[0][0] = msh->coords.data[v1*msh->coords.ncols];
        xs[0][1] = msh->coords.data[v1*msh->coords.ncols+1];
        xs[1][0] = msh->coords.data[v2*msh->coords.ncols];
        xs[1][1] = msh->coords.data[v2*msh->coords.ncols+1];
        xs[2][0] = msh->coords.data[v3*msh->coords.ncols];
        xs[2][1] = msh->coords.data[v3*msh->coords.ncols+1];
        if (inside_tri(xs,ps)){
            stop = true;
            return true;
        } else {
            double AB[2] = {xs[1][0]-xs[0][0], xs[1][1]-xs[0][1]};
            double BC[2] = {xs[2][0]-xs[1][0], xs[2][1]-xs[1][1]};
            double CA[2] = {xs[0][0]-xs[2][0], xs[0][1]-xs[2][1]};
            double AP[2] = {ps[0]-xs[0][0], ps[1]-xs[0][1]};
            double BP[2] = {ps[0]-xs[1][0], ps[1]-xs[1][1]};
            double CP[2] = {ps[0]-xs[2][0], ps[1]-xs[2][1]};
            double N1[2] = {AB[1],-AB[0]};
            double N2[2] = {BC[1],-BC[0]};
            double N3[2] = {CA[1],-CA[0]};
            double S1 = AP[0]*N1[0]+AP[1]*N1[1];
            double S2 = BP[0]*N2[0]+BP[1]*N2[1];
            double S3 = CP[0]*N3[0]+CP[1]*N3[1];
            if ((S1>0)&&(S1>=S2)&&(S1>=S3)){
                hfid = msh->sibhfs.data[(*tri)*msh->sibhfs.ncols];
                if (hfid != 0){
                    *tri = hfid2eid(hfid)-1;
                } else {
                    stop = true;
                    *tri = elids2hfid(*tri+1,1);
                    return false;
                }
            } else if ((S2>0)&&(S2>=S1)&&(S2>=S3)) {
                hfid = msh->sibhfs.data[(*tri)*msh->sibhfs.ncols+1];
                if (hfid != 0){
                    *tri = hfid2eid(hfid)-1;
                } else {
                    stop = true;
                    *tri =  elids2hfid(*tri+1,2);
                    return false;
                }
            } else if ((S3>0)&&(S3>=S1)&&(S3>=S2)){
                hfid = msh->sibhfs.data[(*tri)*msh->sibhfs.ncols+2];
                if (hfid != 0){
                    *tri = hfid2eid(hfid)-1;
                } else {
                    stop = true;
                    *tri = elids2hfid(*tri+1,3);
                    return false;
                }
            } else {   
                printf("error occurd\n");
            }
        }
        iters++;
    }
    if (iters >= 500){
        printf("infinite loop encountered\n");
    }
    *tri = -1;
    return false;
}

void flip_edge(struct Mesh* msh, int eid, int lid){
    int hfid = msh->sibhfs.data[3*eid + lid];
    int oppeid = hfid2eid(hfid)-1;
    int opplid = hfid2lid(hfid)-1;

    int v1 = msh->elems.data[3*eid + ((lid+2)%3)];
    int v2 = msh->elems.data[3*eid + lid];
    int v3 = msh->elems.data[3*eid + ((lid+1)%3)];
    int v4 = msh->elems.data[3*oppeid + ((opplid+2)%3)];

    int tri1[3] = {v1,v2,v4};
    int tri2[3] = {v1,v4,v3};
    int sib1[3] = {msh->sibhfs.data[3*eid + ((lid+2)%3)],msh->sibhfs.data[3*oppeid + (opplid+1)%3],elids2hfid(oppeid+1,1)};
    int sib2[3] = {elids2hfid(eid+1,3),msh->sibhfs.data[oppeid*3 + ((opplid+2)%3)],msh->sibhfs.data[eid*3 + ((lid+1)%3)]};

    msh->elems.data[eid*3] = tri1[0]; msh->elems.data[eid*3+1] = tri1[1]; msh->elems.data[eid*3+2] = tri1[2];
    msh->elems.data[oppeid*3] = tri2[0]; msh->elems.data[oppeid*3+1] = tri2[1]; msh->elems.data[oppeid*3+2] = tri2[2];
    msh->sibhfs.data[eid*3] = sib1[0];
    msh->sibhfs.data[eid*3+1] = sib1[1];
    msh->sibhfs.data[eid*3+2] = sib1[2];
    msh->sibhfs.data[oppeid*3] = sib2[0];
    msh->sibhfs.data[oppeid*3+1] = sib2[1];
    msh->sibhfs.data[oppeid*3+2] = sib2[2];

    bool sib11 = false;
    bool sib12 = false;
    bool sib21 = false;
    bool sib22 = false;
    hfid = msh->sibhfs.data[eid*3];
    if (hfid2eid(hfid) > 0){
        msh->sibhfs.data[3*(hfid2eid(hfid)-1) + (hfid2lid(hfid)-1)] = elids2hfid(eid+1,1);
    }
    hfid = msh->sibhfs.data[eid*3+1];
    if (hfid2eid(hfid) > 0){
        msh->sibhfs.data[3*(hfid2eid(hfid)-1) + (hfid2lid(hfid)-1)] = elids2hfid(eid+1,2);
    }
    hfid = msh->sibhfs.data[oppeid*3+1];
    if (hfid2eid(hfid) > 0){
        msh->sibhfs.data[3*(hfid2eid(hfid)-1) + (hfid2lid(hfid)-1)] = elids2hfid(oppeid+1,2);
    }
    hfid = msh->sibhfs.data[oppeid*3+2];
    if (hfid2eid(hfid) > 0){
        msh->sibhfs.data[3*(hfid2eid(hfid)-1) + (hfid2lid(hfid)-1)] = elids2hfid(oppeid+1,3);
    }
    return;
}

bool inside_tri(const double xs[3][2], const double ps[2]){
    double val1 = (ps[0]-xs[1][0])*(xs[0][1]-xs[1][1]) - (xs[0][0]-xs[1][0])*(ps[1]-xs[1][1]);
    double val2 = (ps[0]-xs[2][0])*(xs[1][1]-xs[2][1]) - (xs[1][0]-xs[2][0])*(ps[1]-xs[2][1]);
    double val3 = (ps[0]-xs[0][0])*(xs[2][1]-xs[0][1]) - (xs[2][0]-xs[0][0])*(ps[1]-xs[0][1]);
    bool has_neg = (val1 <= 0) || (val2<=0) || (val3<=0);
    bool has_pos = (val1 >= 0) || (val2>=0) || (val3>=0);
    return !(has_neg && has_pos);
}

bool inside_circumtri(const double xs[3][2], const double ps[2]){
    double* C = circumcenter(xs);
    double R = (xs[0][0]-C[0])*(xs[0][0]-C[0]) + (xs[0][1] - C[1])*(xs[0][1] - C[1]);
    bool D = ((ps[0]-C[0])*(ps[0]-C[0]) + (ps[1] - C[1])*(ps[1] - C[1])) < R;
    return (D);
}

double* circumcenter(const double xs[3][2]){
    double ax = xs[0][0];
    double ay = xs[0][1];
    double bx = xs[1][0];
    double by = xs[1][1];
    double cx = xs[2][0];
    double cy = xs[2][1];
    double D = 2*(ax*(by-cy) + bx*(cy-ay) + cx*(ay-by));
    double ux = (ax*ax + ay*ay)*(by-cy) + \
        (bx*bx + by*by)*(cy-ay) + \
        (cx*cx + cy*cy)*(ay-by);
    double uy = (ax*ax + ay*ay)*(cx-bx) + \
        (bx*bx + by*by)*(ax-cx) + \
        (cx*cx + cy*cy)*(bx-ax);
    static double C[2];
    C[0] = ux/D; C[1] = uy/D;
    return C;
}

bool Line_cross(double p1x, double p1y, double p2x, double p2y, double p3x, double p3y, double p4x, double p4y){
    bool a1 = (p4y-p1y)*(p3x-p1x) > (p3y-p1y)*(p4x-p1x);
    bool a2 = (p4y-p2y)*(p3x-p2x) > (p3y-p2y)*(p4x-p2x);
    bool b1 = (p3y-p1y)*(p2x-p1x) > (p2y-p1y)*(p3x-p1x);
    bool b2 = (p4y-p1y)*(p2x-p1x) > (p2y-p1y)*(p4x-p1x);

    return (a1 != a2 && b1 != b2);
}

void Mesh_draw(struct Mesh* msh){
    int nv = msh->coords.nrows;
    double x[nv];
    double y[nv];
    for (int i = 0; i<nv; i++){
        x[i] = msh->coords.data[i*msh->coords.ncols];
        y[i] = msh->coords.data[i*msh->coords.ncols+1];
    }

    bool success;
	StringReference *errorMessage;
	RGBABitmapImageReference *imageReference = CreateRGBABitmapImageReference();
	ScatterPlotSeries *series = GetDefaultScatterPlotSeriesSettings();
	series->xs = x;
	series->xsLength = sizeof(x)/sizeof(double);
	series->ys = y;
	series->ysLength = sizeof(y)/sizeof(double);
	series->linearInterpolation = false;
	series->pointType = L"dots";
	series->pointTypeLength = wcslen(series->pointType);
    series->lineThickness = 1;
	series->color = CreateRGBColor(1, 0, 0);

	ScatterPlotSettings *settings = GetDefaultScatterPlotSettings();
	settings->width = 500;
	settings->height = 500;
	settings->autoBoundaries = true;
	settings->autoPadding = true;
	settings->title = L"";
	settings->titleLength = wcslen(settings->title);
	settings->xLabel = L"";
	settings->xLabelLength = wcslen(settings->xLabel);
	settings->yLabel = L"";
	settings->yLabelLength = wcslen(settings->yLabel);
	ScatterPlotSeries *s [] = {series};
	settings->scatterPlotSeries = s;
	settings->scatterPlotSeriesLength = 1;
    settings->showGrid = false;
    settings->autoBoundaries = false;
    settings->autoPadding = false;
    settings->xMax = 0.8;
    settings->xMin = 0.2;
    settings->yMax = 0.8;
    settings->yMin = 0.2;


	errorMessage = (StringReference *)malloc(sizeof(StringReference));
	success = DrawScatterPlotFromSettings(imageReference, settings, errorMessage);

    double x1,x2,y1,y2;
    int ncols = msh->elems.ncols;
    for (int i = 0; i<msh->nelems; i++){
        for (int j = 0; j<ncols; j++){
            x1 = MapXCoordinateBasedOnSettings(msh->coords.data[2*msh->elems.data[ncols*i + j]], settings);
            y1 = MapYCoordinateBasedOnSettings(msh->coords.data[2*msh->elems.data[ncols*i + j]+1], settings);
            x2 = MapXCoordinateBasedOnSettings(msh->coords.data[2*msh->elems.data[ncols*i + (j+1)%ncols]], settings);
            y2 = MapYCoordinateBasedOnSettings(msh->coords.data[2*msh->elems.data[ncols*i + (j+1)%ncols]+1], settings);
            DrawLine(imageReference->image, x1,y1,x2,y2,1, CreateRGBColor(0, 0, 0));
        }
    }


    if(success){
		size_t length;
		double *pngdata = ConvertToPNG(&length, imageReference->image);
		WriteToFile(pngdata, length, "mesh.png");
		DeleteImage(imageReference->image);
	}
    FreeAllocations();

    return;
}

