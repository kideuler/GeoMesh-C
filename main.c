#include "GeoMesh.h"

int main(void){
    srand ( time(NULL) );
    int npoints = 100;

    struct DoubleMatrix xs;
    struct IntMatrix segs;
    struct DoubleMatrix ps;
    struct IntMatrix segs2;
    double h = Ellipse(&segs, &xs, npoints,true);
    int nv = 3*npoints;
    double h2 = Flower(&segs2, &ps, nv, false);
    //load_lakeSuperior(&segs, &xs);
    //load_monalisa(&ps);
    nv = ps.nrows;
    
    clock_t start = clock(), diff;
    //struct Mesh msh = GeoMesh_Delaunay(&xs);
    struct Mesh msh = GeoMesh_ConstrainedDelaunay(&segs,&xs);
    diff = clock() - start;
    int msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Created %d Triangles from initial in %d seconds %d milliseconds\n", msh.nelems, msec/1000, msec%1000);

    double* h_ratios = (double*) malloc(nv*sizeof(double));
    for (int i = 0; i<nv; i++){
        h_ratios[i] = 0.1;
    }
    start = clock();
    msh.haskdTree = false;
    //Mesh_compute_kdTree(&msh, ps, h, h_ratios, 0.1);
    diff = clock() - start;
    msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Created kd-Tree with %d nodes in %d seconds %d milliseconds\n", nv, msec/1000, msec%1000);

    bool is_sib = check_sibhfs(&msh);

    start = clock();
    GeoMesh_DelaunayRefine(&msh,true, h,2);
    diff = clock() - start;
    msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Created %d Triangles from refinement in %d seconds %d milliseconds\n", msh.nelems, msec/1000, msec%1000);

    bool* bdy = Mesh_find_bdy_nodes(&msh);

    start = clock();
    Mesh_smooth2d(&msh, bdy, 500);
    diff = clock() - start;
    msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("finished mesh smoothing in %d seconds %d milliseconds\n", msec/1000, msec%1000);

    GeoMesh_partition(&msh, 1, 2);

    Mesh2vtk(&msh);
    free(bdy);
    return 0;
}