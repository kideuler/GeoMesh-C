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

// Delaunay triangulation driver function
struct Mesh GeoMesh_Delaunay(struct DoubleMatrix *xs){
    int nv = xs->nrows;

    struct Mesh msh;
    struct Mesh *mshbuff = (struct Mesh *) malloc(sizeof(struct Mesh));
    msh.nelems = 0;
    int upper_bound = nv*nv*2;

    // creating temporary buffers
    struct DoubleMatrix *coords = (struct DoubleMatrix *) malloc(sizeof(struct DoubleMatrix));
    coords->nrows = nv+3;
    coords->ncols = 2;
    coords->data = (double* ) malloc((nv+3)*2*sizeof(double));

    mshbuff->elems.nrows = upper_bound;
    mshbuff->elems.ncols = 3;
    mshbuff->elems.data = (int*) malloc((upper_bound)*3*sizeof(int));

    mshbuff->sibhfs.nrows = upper_bound;
    mshbuff->sibhfs.ncols = 3;
    mshbuff->sibhfs.data = (int*) malloc((upper_bound)*3*sizeof(int));


    double ax = DoubleMatrix_min(xs,0); double ay = DoubleMatrix_min(xs,1);
    double bx = DoubleMatrix_max(xs,0); double by = DoubleMatrix_max(xs,1);
    double d = bx-ax > by-bx ? bx-ax : by-ay;

    for (int n = 0; n<nv; n++){
        coords->data[n*coords->ncols] = (xs->data[n*xs->ncols] - ax)/d;
        coords->data[n*coords->ncols+1] = (xs->data[n*xs->ncols+1] - ay)/d;
    }

    // sort points by proximity NOT IMPL.

    // create big triangle
    coords->data[nv*coords->ncols] = -100.0;
    coords->data[nv*coords->ncols+1] = -100.0;
    coords->data[(nv+1)*coords->ncols] = 100.0;
    coords->data[(nv+1)*coords->ncols+1] = -100.0;
    coords->data[(nv+2)*coords->ncols] = 0.0;
    coords->data[(nv+2)*coords->ncols+1] = 100.0;

    mshbuff->elems.data[0] = nv;
    mshbuff->elems.data[1] = nv+1; 
    mshbuff->elems.data[2] = nv+2;
    mshbuff->sibhfs.data[0] = 0; mshbuff->sibhfs.data[1] = 0; mshbuff->sibhfs.data[2] = 0;

    // main loop inserting each node into the triangulation
    int enclosed_tri=-1;
    int vid;
    bool exitl, inside;
    mshbuff->coords = *coords;
    mshbuff->nelems = 1;
    mshbuff->delete_elem = (bool *) malloc(upper_bound*sizeof(bool));
    for (int n = 0; n<nv; n++){
        vid = n;
        enclosed_tri = mshbuff->nelems-1;
        // find enclosing triangle
        double ps[2] = {mshbuff->coords.data[vid*mshbuff->coords.ncols],mshbuff->coords.data[vid*mshbuff->coords.ncols+1]};
        inside = Mesh_find_enclosing_tri(mshbuff, &enclosed_tri, ps);

        if (!inside){ printf("no enclosing triangle found\n");}

        // insert node into triangulation using sloans flip insertion algorithm
        Mesh_flip_insertion(mshbuff, &vid, enclosed_tri);
    }




    msh.coords = *xs;
    msh.nelems = mshbuff->nelems;
    msh.elems = IntMatrix_create(msh.nelems,3);
    msh.sibhfs = IntMatrix_create(msh.nelems,3);
    int k = 0;
    for (int i = 0; i<msh.nelems; i++){
        if (!mshbuff->delete_elem[i]) {
            for (int j = 0; j<3; j++){
                msh.elems.data[msh.elems.ncols*k + j] = mshbuff->elems.data[msh.elems.ncols*i + j];
                msh.sibhfs.data[msh.sibhfs.ncols*k + j] = mshbuff->sibhfs.data[msh.sibhfs.ncols*i + j];
            }
            k++;
        }
    }

    // freeing temporary buffers
    free(coords->data); free(coords); 
    free(mshbuff->elems.data);
    free(mshbuff->sibhfs.data);
    free(mshbuff->delete_elem);
    free(mshbuff);
    return msh;
}

void Mesh_flip_insertion(struct Mesh* msh, int* vid, int tri_start){
    msh->delete_elem[tri_start] = true;
    int hfid,eid,lid;

    int tri[3] = {msh->elems.data[tri_start*msh->elems.ncols],msh->elems.data[tri_start*msh->elems.ncols+1],msh->elems.data[tri_start*msh->elems.ncols+2]};
    int sib[3] = {msh->sibhfs.data[tri_start*msh->sibhfs.ncols],msh->sibhfs.data[tri_start*msh->sibhfs.ncols+1],msh->sibhfs.data[tri_start*msh->sibhfs.ncols+2]};

    int* stack = (int*) malloc(msh->nelems*sizeof(int));
    int stack_size = 0;
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
            stack[stack_size] = elids2hfid(eids[i]+1, 2);
            stack_size++;
        }
    }

    msh->nelems+=3;
    double xs[3][2];
    int oppeid, opplid;
    while (stack_size > -1){
        hfid = stack[stack_size];
        stack_size--;
        eid = hfid2eid(hfid)-1;
        lid = hfid2lid(hfid)-1;
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
    while (!stop){
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

            }
        }
        iters++;
    }
    *tri = -1;
    return false;
}

bool inside_tri(const double xs[3][2], const double ps[2]){
    double val1 = (ps[0]-xs[1][0])*(xs[0][1]-xs[1][1]) - (xs[0][0]-xs[1][0])*(ps[1]-xs[1][1]);
    double val2 = (ps[0]-xs[2][0])*(xs[1][1]-xs[2][1]) - (xs[1][0]-xs[2][0])*(ps[1]-xs[2][1]);
    double val3 = (ps[0]-xs[0][0])*(xs[2][1]-xs[0][1]) - (xs[2][0]-xs[0][0])*(ps[1]-xs[0][1]);
    bool has_neg = (val1 <= 0) || (val2<=0) || (val3<=0);
    bool has_pos = (val1 >= 0) || (val2>=0) || (val3>=0);
    return !(has_neg && has_pos);
}