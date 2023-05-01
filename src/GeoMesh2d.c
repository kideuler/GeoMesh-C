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

// Use Delaunay mesh to refine mesh
void GeoMesh_DelaunayRefine(struct Mesh *msh, bool use_edgelengh, double h_target, int point_algorithm){
    int nv = msh->coords.nrows;
    int nelems_initial = msh->nelems;
    int n,i,tri,hfid,eid,lid,upper_bound,e;
    bool exitl;
    double alpha=1.0;
    double theta,shape;
    double ps[3][2];

    // finding total area to estimate size bounds
    double total_area = 0.0;
    double area;
    for (n=0; n<nelems_initial; n++){
        ps[0][0] = msh->coords.data[2*msh->elems.data[3*n]];
        ps[0][1] = msh->coords.data[2*msh->elems.data[3*n]+1];
        ps[1][0] = msh->coords.data[2*msh->elems.data[3*n+1]];
        ps[1][1] = msh->coords.data[2*msh->elems.data[3*n+1]+1];
        ps[2][0] = msh->coords.data[2*msh->elems.data[3*n+2]];
        ps[2][1] = msh->coords.data[2*msh->elems.data[3*n+2]+1];
        area = area_tri(ps);
        total_area += area;
    }

    double radius_target;
    if (use_edgelengh){
        double area_single = pow(sqrt(3)*h_target/3,2)/2;
        radius_target = (sqrt(3)/3)*h_target;
        upper_bound = (int) 2*total_area/area_single;
        upper_bound = (upper_bound>10*msh->nelems)? upper_bound : 10*msh->nelems;
    } else {
        upper_bound = 10*msh->nelems;
    }
    DoubleMatrix_resize(&msh->coords,upper_bound);
    IntMatrix_resize(&msh->elems,upper_bound);
    IntMatrix_resize(&msh->sibhfs,upper_bound);
    free(msh->delete_elem);
    msh->delete_elem = (bool*) malloc(upper_bound*sizeof(bool));
    memset(msh->delete_elem, false, upper_bound*sizeof(bool));
    bool* on_bdy = (bool*) malloc(upper_bound*sizeof(bool));
    memset(on_bdy, false, upper_bound*sizeof(bool));
    for(int i = 0; i<msh->nelems; i++){
        on_bdy[i] = msh->on_boundary[i];
    }
    free(msh->on_boundary);
    msh->on_boundary = on_bdy;
    double *C = (double*) malloc(2*sizeof(double)); // buffer for added center
    int* order = (int*) malloc(upper_bound*sizeof(int));
    for(n=0;n<upper_bound;n++){order[n]=n;}

    // main loop
    n = 0;
    bool inside_domain,stop;
    while (n<msh->nelems){
        e = order[n];
        if (!msh->delete_elem[e]){
            ps[0][0] = msh->coords.data[2*msh->elems.data[3*e]];
            ps[0][1] = msh->coords.data[2*msh->elems.data[3*e]+1];
            ps[1][0] = msh->coords.data[2*msh->elems.data[3*e+1]];
            ps[1][1] = msh->coords.data[2*msh->elems.data[3*e+1]+1];
            ps[2][0] = msh->coords.data[2*msh->elems.data[3*e+2]];
            ps[2][1] = msh->coords.data[2*msh->elems.data[3*e+2]+1];
            if (point_algorithm == 1) {circumcenter(ps, C);}
            if (point_algorithm == 2) {off_circumcenter(ps, sqrt(2.0), C);}
            C[0] += 1e-4*radius_target*drand(-1.0,1.0);
            C[1] += 1e-4*radius_target*drand(-1.0,1.0);
            
            shape = eval_trishape(ps);
            if (use_edgelengh){alpha = eval_alpha(ps,radius_target);}
                
            if (alpha > 1.2 || shape > 2){
                // add coordinate
                msh->coords.data[2*nv] = C[0];
                msh->coords.data[2*nv+1] = C[1];
                tri = e;
                inside_domain = Mesh_find_enclosing_tri(msh,&tri,C);
                if (!inside_domain){
                    if (tri == -1){
                        printf("find triangle location failed: deleting point\n");
                        nv--;
                    } else { // place point on segment midpoint
                        if (inside_diametral(msh,tri, C)){
                            Mesh_flip_insertion_segment(msh,nv, tri);
                        } else {
                            if (msh->on_boundary[e]){
                                Mesh_flip_insertion_segment(msh,nv, tri);
                            } else{
                                if (n == msh->nelems-1 || e == msh->nelems-1){
                                    Mesh_flip_insertion_segment(msh,nv, tri);
                                } else {
                                    nv--;
                                    order[n] = msh->nelems-1;
                                    order[msh->nelems-1] = e;
                                    n--;
                                }
                            }
                        }
                    }
                } else { // check if point is in diametral circle of nearby elems
                    stop = false;
                    if (msh->on_boundary[tri]){
                        hfid = (msh->sibhfs.data[3*tri]==0) ? elids2hfid(tri+1,1) : \
                        ((msh->sibhfs.data[3*tri+1]==0) ? elids2hfid(tri+1,2) : elids2hfid(tri+1,3));
                        if (inside_diametral(msh,hfid, C)){
                            Mesh_flip_insertion_segment(msh,nv, hfid);
                            stop = true;
                        }
                    } else if(msh->sibhfs.data[3*tri] > 0) {
                        if(msh->on_boundary[hfid2eid(msh->sibhfs.data[3*tri])-1]){
                        eid = hfid2eid(msh->sibhfs.data[3*tri])-1;
                        hfid = (msh->sibhfs.data[3*eid]==0) ? elids2hfid(eid+1,1) : \
                        ((msh->sibhfs.data[3*eid+1]==0) ? elids2hfid(eid+1,2) : elids2hfid(eid+1,3));
                        if (inside_diametral(msh,hfid, C)){
                            Mesh_flip_insertion_segment(msh,nv, hfid);
                            stop = true;
                        }
                        }
                    } else if(msh->sibhfs.data[3*tri+1] > 0) {
                        if(msh->on_boundary[hfid2eid(msh->sibhfs.data[3*tri+1])-1]){
                        eid = hfid2eid(msh->sibhfs.data[3*tri+1])-1;
                        hfid = (msh->sibhfs.data[3*eid]==0) ? elids2hfid(eid+1,1) : \
                        ((msh->sibhfs.data[3*eid+1]==0) ? elids2hfid(eid+1,2) : elids2hfid(eid+1,3));
                        if (inside_diametral(msh,hfid, C)){
                            Mesh_flip_insertion_segment(msh,nv, hfid);
                            stop = true;
                        }
                        }
                    } else if(msh->sibhfs.data[3*tri+2] > 0) {
                        if(msh->on_boundary[hfid2eid(msh->sibhfs.data[3*tri+2])-1]){
                        eid = hfid2eid(msh->sibhfs.data[3*tri+2])-1;
                        hfid = (msh->sibhfs.data[3*eid]==0) ? elids2hfid(eid+1,1) : \
                        ((msh->sibhfs.data[3*eid+1]==0) ? elids2hfid(eid+1,2) : elids2hfid(eid+1,3));
                        if (inside_diametral(msh,hfid, C)){
                            Mesh_flip_insertion_segment(msh,nv, hfid);
                            stop = true;
                        }
                        }
                    }

                    if (!stop){
                        Mesh_flip_insertion(msh,&nv,tri);
                    }
                }
                nv++;
            }

        }

        if (msh->nelems > upper_bound-3){
            printf("max buffer size reached, resizing\n");
            upper_bound = 2*upper_bound;
            Mesh_resize(msh, upper_bound);
            free(order);
            order = (int*) malloc(upper_bound*sizeof(int));
            for(i=0;i<upper_bound;i++){order[i]=i;}
            printf("%d %d\n",n,msh->nelems);
        }
        n++;
    }

    if (msh->nelems >= upper_bound-3){
        printf("max buffer size reached\n");
    }
    free(order); free(C);
    Mesh_deleteElems(msh);
    DoubleMatrix_resize(&msh->coords, nv);
}

// Constrained Delaunay triangulation driver function
struct Mesh GeoMesh_ConstrainedDelaunay(struct IntMatrix *segments, struct DoubleMatrix *xs){
    struct Mesh msh = GeoMesh_Delaunay(xs);
    int maxne=20;
    Mesh_compute_OneringElements(&msh,maxne);
    maxne = msh.stncl.maxne;

    int* facets = (int*) malloc(2*msh.nelems*sizeof(int));
    msh.bwork = (bool*) malloc(msh.nelems*sizeof(bool));
    memset(msh.bwork, false, msh.nelems*sizeof(bool));
    msh.delete_elem = (bool*) malloc(msh.nelems*sizeof(bool));
    memset(msh.delete_elem, false, msh.nelems*sizeof(bool));

    int v1,v2,hvid,kk,eid,nid,nf,hfid,lid,oppeid,opplid;
    bool inray = false;
    bool connected;
    bool convex = false;
    int iter = 0;
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
                        msh.bwork[hfid2eid(hfid)-1] = true;
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
                lid = (nid+1)%3;
                hfid = msh.sibhfs.data[3*eid+lid];
                while (!convex_quad(&msh, eid, lid)){ // check if convex quadrilateral
                    oppeid = hfid2eid(hfid)-1;
                    opplid = hfid2lid(hfid)-1;
                    if (Line_cross(msh.coords.data[2*msh.elems.data[3*oppeid + (opplid+1)%3]],msh.coords.data[2*msh.elems.data[3*oppeid + (opplid+1)%3]+1], \
                    msh.coords.data[2*msh.elems.data[3*oppeid + (opplid+2)%3]],msh.coords.data[2*msh.elems.data[3*oppeid + (opplid+2)%3]+1], \
                    msh.coords.data[2*v1],msh.coords.data[2*v1+1], msh.coords.data[2*v2],msh.coords.data[2*v2+1])){
                        eid = oppeid;
                        lid = (opplid+1)%3;
                    } else if (Line_cross(msh.coords.data[2*msh.elems.data[3*oppeid + (opplid+2)%3]],msh.coords.data[2*msh.elems.data[3*oppeid + (opplid+2)%3]+1], \
                    msh.coords.data[2*msh.elems.data[3*oppeid + opplid]],msh.coords.data[2*msh.elems.data[3*oppeid + opplid]+1], \
                    msh.coords.data[2*v1],msh.coords.data[2*v1+1], msh.coords.data[2*v2],msh.coords.data[2*v2+1])){
                        eid = oppeid;
                        lid = (opplid+2)%3;
                    } else {
                        printf("no line crossings found\n");
                    }
                }
                flip_edge(&msh, eid, lid);
                Mesh_compute_OneringElements(&msh, maxne);
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
    msh.on_boundary = (bool* ) malloc(upper_bound*sizeof(bool));
    memset(msh.on_boundary, false, upper_bound*sizeof(bool));
    msh.stack = (int*) malloc(upper_bound*sizeof(int));

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
    coords->data[nv*2] = -100.0;
    coords->data[nv*2+1] = -100.0;
    coords->data[(nv+1)*2] = 100.0;
    coords->data[(nv+1)*2+1] = -100.0;
    coords->data[(nv+2)*2] = 0.0;
    coords->data[(nv+2)*2+1] = 100.0;

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
    memset(msh.delete_elem, false, upper_bound*sizeof(bool));
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
    msh->hasStencil = true;
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
    memset(idx_rev,-1,msh->nelems*sizeof(int));

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
    msh->nelems = nelems;

    for (int i = 0; i<msh->nelems; i++){
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
        if (nside == 3) {
            msh->on_boundary[i] = false;
        } else {
            msh->on_boundary[i] = true;
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
    int hfid,eid,lid;
    msh->delete_elem[tri_start] = true;

    int tri[3] = {msh->elems.data[3*tri_start],msh->elems.data[3*tri_start+1],msh->elems.data[3*tri_start+2]};
    int sib[3] = {msh->sibhfs.data[3*tri_start],msh->sibhfs.data[3*tri_start+1],msh->sibhfs.data[3*tri_start+2]};

    int stack_size = -1;
    int eids[3] = {msh->nelems, msh->nelems+1, msh->nelems+2};
    // splitting triangles and adding them to the stack
    for (int i = 0; i<3; i++){
        msh->elems.data[3*eids[i]] = *vid;
        msh->elems.data[3*eids[i]+1] = tri[i];
        msh->elems.data[3*eids[i]+2] = tri[(i+1)%3];
        hfid = sib[i];
        msh->sibhfs.data[3*eids[i]] = elids2hfid(eids[(i+2)%3]+1 ,3);
        msh->sibhfs.data[3*eids[i]+1] = hfid;
        msh->sibhfs.data[3*eids[i]+2] = elids2hfid(eids[(i+1)%3]+1,1);
        if (hfid2eid(hfid) > 0){
            msh->sibhfs.data[3*(hfid2eid(hfid)-1) + hfid2lid(hfid)-1] = elids2hfid(eids[i]+1, 2);
            msh->stack[stack_size+1] = elids2hfid(eids[i]+1, 2);
            stack_size++;
            msh->on_boundary[eids[i]] = false;
        } else {
            msh->on_boundary[eids[i]] = true;
        }
    }

    msh->nelems+=3;
    double xs[3][2], ps[2];
    int oppeid, opplid;
    while (stack_size > -1){
        hfid = msh->stack[stack_size];
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
                msh->stack[stack_size+1] = elids2hfid(oppeid+1,2);
                stack_size++;
            }
            hfid = msh->sibhfs.data[3*eid+1];
            if (hfid2eid(hfid) > 0){
                msh->stack[stack_size+1] = elids2hfid(eid+1,2);
                stack_size++;
            }
            
        }
        
    }

    return;
}

void Mesh_flip_insertion_segment(struct Mesh* msh, int vid, int hfid){
    int stack_size = -1;
    int eid = hfid2eid(hfid)-1;
    int lid = hfid2lid(hfid)-1;

    msh->coords.data[2*vid] = (msh->coords.data[2*msh->elems.data[3*eid+lid]]+ \
    msh->coords.data[2*msh->elems.data[3*eid+(lid+1)%3]])/2.0;
    msh->coords.data[2*vid+1] = (msh->coords.data[2*msh->elems.data[3*eid+lid]+1]+ \
    msh->coords.data[2*msh->elems.data[3*eid+(lid+1)%3]+1])/2.0;

    msh->delete_elem[eid] = true;
    msh->elems.data[3*msh->nelems] = vid;
    msh->elems.data[3*msh->nelems+1] = msh->elems.data[3*eid+(lid+2)%3];
    msh->elems.data[3*msh->nelems+2] = msh->elems.data[3*eid+lid];
    msh->elems.data[3*(msh->nelems+1)] = vid;
    msh->elems.data[3*(msh->nelems+1)+1] = msh->elems.data[3*eid+(lid+1)%3];
    msh->elems.data[3*(msh->nelems+1)+2] = msh->elems.data[3*eid+(lid+2)%3];

    msh->sibhfs.data[3*msh->nelems] = elids2hfid(msh->nelems+2,3);
    msh->sibhfs.data[3*msh->nelems+1] = msh->sibhfs.data[3*eid+(lid+2)%3];
    msh->sibhfs.data[3*msh->nelems+2] = 0;
    msh->sibhfs.data[3*(msh->nelems+1)] = 0;
    msh->sibhfs.data[3*(msh->nelems+1)+1] = msh->sibhfs.data[3*eid+(lid+1)%3];
    msh->sibhfs.data[3*(msh->nelems+1)+2] = elids2hfid(msh->nelems+1,1);
    if (msh->sibhfs.data[3*eid+(lid+1)%3] != 0){
        msh->sibhfs.data[3*(hfid2eid(msh->sibhfs.data[3*eid+(lid+1)%3])-1)+\
        hfid2lid(msh->sibhfs.data[3*eid+(lid+1)%3])-1] = elids2hfid(msh->nelems+2, 2);
        msh->stack[stack_size+1] = elids2hfid(msh->nelems+2,2);
        stack_size++;
    }
    if (msh->sibhfs.data[3*eid+(lid+2)%3] != 0){
        msh->sibhfs.data[3*(hfid2eid(msh->sibhfs.data[3*eid+(lid+2)%3])-1)+\
        hfid2lid(msh->sibhfs.data[3*eid+(lid+2)%3])-1] = elids2hfid(msh->nelems+1, 2);
        msh->stack[stack_size+1] = elids2hfid(msh->nelems+1,2);
        stack_size++;
    }
    msh->on_boundary[msh->nelems] = true;
    msh->on_boundary[msh->nelems+1] = true;

    msh->nelems+=2;

    double xs[3][2],ps[2];
    int oppeid, opplid;
    while (stack_size > -1){
        hfid = msh->stack[stack_size];
        stack_size--;
        eid = hfid2eid(hfid)-1;
        lid = hfid2lid(hfid)-1;
        oppeid = hfid2eid(msh->sibhfs.data[3*eid+lid])-1;
        opplid = hfid2lid(msh->sibhfs.data[3*eid+lid])-1;
        xs[0][0] = msh->coords.data[2*msh->elems.data[3*oppeid]];
        xs[0][1] = msh->coords.data[2*msh->elems.data[3*oppeid]+1];
        xs[1][0] = msh->coords.data[2*msh->elems.data[3*oppeid+1]];
        xs[1][1] = msh->coords.data[2*msh->elems.data[3*oppeid+1]+1];
        xs[2][0] = msh->coords.data[2*msh->elems.data[3*oppeid+2]];
        xs[2][1] = msh->coords.data[2*msh->elems.data[3*oppeid+2]+1];
        ps[0] = msh->coords.data[2*msh->elems.data[3*eid]];
        ps[1] = msh->coords.data[2*msh->elems.data[3*eid]+1];

        if (inside_circumtri(xs,ps)){
            flip_edge(msh, eid,lid);

            hfid = msh->sibhfs.data[3*oppeid+1];
            if (hfid2eid(hfid) > 0){
                msh->stack[stack_size+1] = elids2hfid(oppeid+1,2);
                stack_size++;
            }
            hfid = msh->sibhfs.data[3*eid+1];
            if (hfid2eid(hfid) > 0){
                msh->stack[stack_size+1] = elids2hfid(eid+1,2);
                stack_size++;
            }
        }
    }
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
                //assert(false);
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
    } else {
        sib11 = true;
    }
    hfid = msh->sibhfs.data[eid*3+1];
    if (hfid2eid(hfid) > 0){
        msh->sibhfs.data[3*(hfid2eid(hfid)-1) + (hfid2lid(hfid)-1)] = elids2hfid(eid+1,2);
    } else {
        sib12 = true;
    }
    hfid = msh->sibhfs.data[oppeid*3+1];
    if (hfid2eid(hfid) > 0){
        msh->sibhfs.data[3*(hfid2eid(hfid)-1) + (hfid2lid(hfid)-1)] = elids2hfid(oppeid+1,2);
    } else {
        sib21 = true;
    }
    hfid = msh->sibhfs.data[oppeid*3+2];
    if (hfid2eid(hfid) > 0){
        msh->sibhfs.data[3*(hfid2eid(hfid)-1) + (hfid2lid(hfid)-1)] = elids2hfid(oppeid+1,3);
    } else {
        sib22 = true;
    }
    msh->on_boundary[eid] = sib11 || sib12;
    msh->on_boundary[oppeid] = sib21 || sib22;

    if (false){
        // fix stencil for v1
        int nhvid = msh->stncl.nhvids[v1];
        for (int i = 0; i<nhvid; i++){
            if (hfid2eid(msh->stncl.hvids.data[v1*msh->stncl.maxne + i])-1 == eid){
                msh->stncl.hvids.data[v1*msh->stncl.maxne + i] = elids2hfid(eid+1,1);
                msh->stncl.hvids.data[v1*msh->stncl.maxne + nhvid] = elids2hfid(oppeid+1,1);
                msh->stncl.nhvids[v1]++;
                break;
            }
        }

        // fix stencil for v2
        int place,kk;
        nhvid = msh->stncl.nhvids[v2];
        for (int i = 0; i<nhvid; i++){
            if (hfid2eid(msh->stncl.hvids.data[v2*msh->stncl.maxne + i])-1 == eid){
                msh->stncl.hvids.data[v2*msh->stncl.maxne + i] = elids2hfid(eid+1,2);
                break;
            }
            if (hfid2eid(msh->stncl.hvids.data[v2*msh->stncl.maxne + i])-1 == oppeid){
                place = i;
                break;
            }
        }
        kk=0;
        for (int i = 0; i<nhvid; i++){
            if (i != place){
                msh->stncl.hvids.data[v2*msh->stncl.maxne + kk] = msh->stncl.hvids.data[v2*msh->stncl.maxne + i];
                kk++;
            }
        }
        msh->stncl.nhvids[v1]--;

        // fix stencil for v3
        nhvid = msh->stncl.nhvids[v3];
        for (int i = 0; i<nhvid; i++){
            if (hfid2eid(msh->stncl.hvids.data[v3*msh->stncl.maxne + i])-1 == oppeid){
                msh->stncl.hvids.data[v3*msh->stncl.maxne + i] = elids2hfid(oppeid+1,3);
                break;
            }
            if (hfid2eid(msh->stncl.hvids.data[v3*msh->stncl.maxne + i])-1 == eid){
                place = i;
                break;
            }
        }
        kk=0;
        for (int i = 0; i<nhvid; i++){
            if (i != place){
                msh->stncl.hvids.data[v3*msh->stncl.maxne + kk] = msh->stncl.hvids.data[v3*msh->stncl.maxne + i];
                kk++;
            }
        }
        msh->stncl.nhvids[v1]--;

        // fix stencil for v4
        nhvid = msh->stncl.nhvids[v4];
        for (int i = 0; i<nhvid; i++){
            if (hfid2eid(msh->stncl.hvids.data[v4*msh->stncl.maxne + i])-1 == oppeid){
                msh->stncl.hvids.data[v4*msh->stncl.maxne + i] = elids2hfid(oppeid+1,2);
                msh->stncl.hvids.data[v4*msh->stncl.maxne + nhvid] = elids2hfid(eid+1,3);
                msh->stncl.nhvids[v1]++;
                break;
            }
        }
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
    
    double* C = (double*) malloc(2*sizeof(double));
    circumcenter(xs, C);
    double R = (xs[0][0]-C[0])*(xs[0][0]-C[0]) + (xs[0][1] - C[1])*(xs[0][1] - C[1]);
    bool D = ((ps[0]-C[0])*(ps[0]-C[0]) + (ps[1] - C[1])*(ps[1] - C[1])) < R;
    free(C);
    return (D);
}

bool inside_diametral(struct Mesh* msh, int hfid, double ps[2]){
    int eid = hfid2eid(hfid) - 1;
    int lid = hfid2lid(hfid) - 1;
    double v1[2] = {msh->coords.data[2*msh->elems.data[3*eid+lid]],msh->coords.data[2*msh->elems.data[3*eid+lid]+1]};
    double v2[2] = {msh->coords.data[2*msh->elems.data[3*eid+(lid+1)%3]],msh->coords.data[2*msh->elems.data[3*eid+(lid+1)%3]+1]};
    double r = sqrt(pow(v2[0]-v1[0],2) + pow(v2[1]-v1[1],2))/2;
    double p[2] = {(v1[0]+v2[0])/2, (v1[1]+v2[1])/2}; 
    double dist = sqrt(pow(ps[0]-p[0],2) + pow(ps[1]-p[1],2));
    return dist < r;
}

void circumcenter(const double xs[3][2], double* C){
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
    C[0] = ux/D; C[1] = uy/D;
    return;
}

void off_circumcenter(const double xs[3][2], double beta, double* C){

    double* c1 = (double*) malloc(2*sizeof(double));
    circumcenter(xs, c1);
    double ps[3][2];
    double m[2];
    double distpq = 1.0e6;
    double temp=0.0;
    for (int i = 0; i<3; i++){
        temp = sqrt(pow(xs[(i+1)%3][0]-xs[i][0],2)+pow(xs[(i+1)%3][1]-xs[i][1],2));
        if (temp < distpq){
            distpq = temp;
            m[0] = (xs[(i+1)%3][0]+xs[i][0])/2;
            m[1] = (xs[(i+1)%3][1]+xs[i][1])/2;
            ps[0][0] = xs[i][0];
            ps[0][1] = xs[i][1];
            ps[1][0] = xs[(i+1)%3][0];
            ps[1][1] = xs[(i+1)%3][1];
        }
    }
    ps[2][0] = c1[0];
    ps[2][1] = c1[1];
    double* c2 = (double*) malloc(2*sizeof(double));
    circumcenter(ps,c2);
    double distc1c2 = sqrt(pow(c1[0]-c2[0],2.0) + pow(c1[1]-c2[1],2.0));
    if (distc1c2 <= beta*distpq){
        C[0] = c1[0];
        C[1] = c1[1];
    } else{
        C[0] = c2[0] + 0.95*beta*distpq*(c2[0]-m[0])/sqrt(pow(c2[0]-m[0],2) + pow(c2[1]-m[1],2));
        C[1] = c2[1] + 0.95*beta*distpq*(c2[1]-m[1])/sqrt(pow(c2[0]-m[0],2) + pow(c2[1]-m[1],2));
    }
    free(c1); free(c2);
    return;
}

bool* Mesh_find_bdy_nodes(struct Mesh* msh){
    int nv = msh->coords.nrows;
    bool* bdy = (bool*) malloc(nv*sizeof(bool));
    memset(bdy, false, nv*sizeof(bool));

    for (int i = 0; i<msh->nelems; i++){
        for (int j = 0; j<3; j++){
            if (msh->sibhfs.data[3*i+j] == 0){
                bdy[msh->elems.data[3*i+j]] = true;
                bdy[msh->elems.data[3*i+(j+1)%3]] = true;
            }
        }
    }
    return bdy;
}

double area_tri(double xs[3][2]){
    double N = (xs[1][0]-xs[0][0])*(xs[2][1]-xs[0][1]) - (xs[1][1]-xs[0][1])*(xs[2][0]-xs[0][0]);
    return N/2;
}

bool Line_cross(double p1x, double p1y, double p2x, double p2y, double p3x, double p3y, double p4x, double p4y){
    bool a1 = (p4y-p1y)*(p3x-p1x) > (p3y-p1y)*(p4x-p1x);
    bool a2 = (p4y-p2y)*(p3x-p2x) > (p3y-p2y)*(p4x-p2x);
    bool b1 = (p3y-p1y)*(p2x-p1x) > (p2y-p1y)*(p3x-p1x);
    bool b2 = (p4y-p1y)*(p2x-p1x) > (p2y-p1y)*(p4x-p1x);

    return (a1 != a2 && b1 != b2);
}

bool convex_quad(struct Mesh* msh, int eid, int lid){
    int hfid = msh->sibhfs.data[3*eid+lid];
    int oppeid = hfid2eid(hfid)-1;
    int opplid = hfid2lid(hfid)-1;

    int v1 = msh->elems.data[3*eid+(lid+2)%3];
    int v2 = msh->elems.data[3*eid+lid];
    int v4 = msh->elems.data[3*eid+(lid+1)%3];
    int v3 = msh->elems.data[3*oppeid+(opplid+2)%3];

    double sides[4][2] = {{msh->coords.data[2*v2]-msh->coords.data[2*v1],msh->coords.data[2*v2+1]-msh->coords.data[2*v1+1]},
    {msh->coords.data[2*v3]-msh->coords.data[2*v2],msh->coords.data[2*v3+1]-msh->coords.data[2*v2+1]},
    {msh->coords.data[2*v4]-msh->coords.data[2*v3],msh->coords.data[2*v4+1]-msh->coords.data[2*v3+1]},
    {msh->coords.data[2*v1]-msh->coords.data[2*v4],msh->coords.data[2*v1+1]-msh->coords.data[2*v4+1]}};

    double cp;
    bool sign = false;
    for (int i = 0; i<4; i++){
        cp = sides[i][0]*sides[(i+1)%4][1] - sides[i][1]*sides[(i+1)%4][0];
        if (i==0){
            sign = (cp>0);
        } else if (sign != (cp>0)){
            return false;
        }
    }

    return true;
}

double eval_alpha(const double xs[3][2],double radius_target){
    double a = sqrt(pow(xs[1][0]-xs[0][0],2)+pow(xs[1][1]-xs[0][1],2));
    double b = sqrt(pow(xs[2][0]-xs[1][0],2)+pow(xs[2][1]-xs[1][1],2));
    double c = sqrt(pow(xs[2][0]-xs[0][0],2)+pow(xs[2][1]-xs[0][1],2));
    double s = 0.5*(a+b+c);
    double A = sqrt(s*(s-a)*(s-b)*(s-c));
    double r = a*b*c/(4*A);
    return r/radius_target;
}

double eval_trishape(const double xs[3][2]){
    double a = sqrt(pow(xs[1][0]-xs[0][0],2) + pow(xs[1][1]-xs[0][1],2));
    double b = sqrt(pow(xs[2][0]-xs[1][0],2) + pow(xs[2][1]-xs[1][1],2));
    double c = sqrt(pow(xs[0][0]-xs[2][0],2) + pow(xs[0][1]-xs[2][1],2));
    double sa = ((xs[1][0]-xs[0][0])*(xs[2][0]-xs[0][0]) + (xs[1][1]-xs[0][1])*(xs[2][1]-xs[0][1]))/(a*c);
    double sb = ((xs[0][0]-xs[1][0])*(xs[2][0]-xs[1][0]) + (xs[0][1]-xs[1][1])*(xs[2][1]-xs[1][1]))/(a*b);
    double sc = ((xs[0][0]-xs[2][0])*(xs[1][0]-xs[2][0]) + (xs[0][1]-xs[2][1])*(xs[1][1]-xs[2][1]))/(b*c);
    sa = sqrt(1-sa*sa); sb = sqrt(1-sb*sb); sc = sqrt(1-sc*sc);
    return (sa+sb+sc)/(4*sa*sb*sc);
}

void Mesh_draw(struct Mesh* msh){
    double xmin = DoubleMatrix_min(&msh->coords,0);
    double ymin = DoubleMatrix_min(&msh->coords,1);
    double xmax = DoubleMatrix_max(&msh->coords,0);
    double ymax = DoubleMatrix_max(&msh->coords,1);
    double xmid = (xmin+xmax)/2; double ymid = (ymin+ymax)/2;
    double d = (xmax-xmin) > (ymax-ymin) ? (xmax-xmin) : (ymax-ymin);
    xmin = xmid-d/2;
    xmax = xmid+d/2;
    ymin = ymid-d/2;
    ymax = ymid+d/2;

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
	settings->width = 3000;
	settings->height = 3000;
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
    settings->xMax = xmax;
    settings->xMin = xmin;
    settings->yMax = ymax;
    settings->yMin = ymin;


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

bool check_sibhfs(struct Mesh* msh){
    const int edges[3][2] = {{0,1},{1,2},{2,0}};
    int nelems = msh->nelems;
    int hfid,eid,lid;
    bool check = true;
    for (int i = 0; i<nelems; i++){
        for (int j = 0; j<3; j++){
            hfid = msh->sibhfs.data[3*i+j];
            if (hfid != 0){
            eid = hfid2eid(hfid);
            lid = hfid2lid(hfid);
            if (eid > nelems || lid > 3 || eid < 0 || msh->delete_elem[eid-1]){
                printf("sibhfs is wrong at elem %d lid %d oppeid %d opplid %d\n",i,j,eid-1,lid-1);
                check = false;
            } 
            if (msh->elems.data[3*i+edges[j][0]] != msh->elems.data[3*(eid-1)+edges[lid-1][1]] || msh->elems.data[3*i+edges[j][1]] != msh->elems.data[3*(eid-1)+edges[lid-1][0]]){
                printf("sides dont match, wrong at elem %d, lid %d, oppeid %d, opplid %d, sibhfs %d\n",i,j,eid-1,lid-1,hfid);
                check = false;
            }   
            }
        }
    }
    return check;
} 

bool check_jacobians(struct Mesh* msh){
    static double dphi[3][2] = {{-1.0,-1.0},{1.0,0.0},{0.0,1.0}};
    double xs[3][2], J[2][2], detJ;
    bool check = true;

    for(int i =0; i<msh->nelems; i++){
        xs[0][0] = msh->coords.data[2*msh->elems.data[3*i]];
        xs[0][1] = msh->coords.data[2*msh->elems.data[3*i]+1];
        xs[1][0] = msh->coords.data[2*msh->elems.data[3*i+1]];
        xs[1][1] = msh->coords.data[2*msh->elems.data[3*i+1]+1];
        xs[2][0] = msh->coords.data[2*msh->elems.data[3*i+2]];
        xs[2][1] = msh->coords.data[2*msh->elems.data[3*i+2]+1];

        for (int ii = 0; ii<2; ii++){
            for (int jj = 0; jj<2; jj++){
                J[ii][jj] = 0.0;
                for (int kk = 0; kk<3; kk++){
                    J[ii][jj] += xs[kk][ii]*dphi[kk][jj];
                }
            }
        }

        detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];
        if (area_tri(xs)<0 || detJ < 0){
            printf("negative jacobian at eid: %d area: %f\n",i,area_tri(xs));
            check = false;
        }
    }
    return check;
}

