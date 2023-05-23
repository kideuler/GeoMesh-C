#include "GeoMesh.h"

int main(void){
    srand ( time(NULL) );
    /*
    int npoints = 50;
    struct DoubleMatrix xs;
    struct IntMatrix segs;
    struct DoubleMatrix ps;
    struct IntMatrix segs2;
    double h = Circle(&segs, &xs, npoints,false);
    int nv = 3*npoints;
    double h2 = Flower(&segs2, &ps, nv, false);
    //load_lakeSuperior(&segs, &xs);
    //load_monalisa(&ps);
    nv = ps.nrows;
    
    clock_t start = clock(), diff;
    struct Mesh msh = GeoMesh_Delaunay(&xs,2);
    //struct Mesh msh = GeoMesh_ConstrainedDelaunay(&segs,&xs);
    diff = clock() - start;
    int msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Created %d Triangles from initial in %d seconds %d milliseconds\n", msh.nelems, msec/1000, msec%1000);

    double* h_ratios = (double*) malloc(nv*sizeof(double));
    for (int i = 0; i<nv; i++){
        h_ratios[i] = 0.05;
    }
    start = clock();
    msh.haskdTree = false;
    Mesh_compute_kdTree(&msh, ps, h, h_ratios, 0.1);
    diff = clock() - start;
    msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Created kd-Tree with %d nodes in %d seconds %d milliseconds\n", nv, msec/1000, msec%1000);

    bool is_sib = check_sibhfs(&msh);

    start = clock();
    GeoMesh_DelaunayRefine(&msh,true, h,2);
    diff = clock() - start;
    msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Created %d Triangles from refinement in %d seconds %d milliseconds\n", msh.nelems, msec/1000, msec%1000);
    int nparts = 4;
    int type = 2;
    Mesh_Graphinit(&msh, type);
    

    start = clock();
    GeoMesh_partition(&msh, type, nparts);
    diff = clock() - start;
    msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("finished mesh patition into %d pieces in %d seconds %d milliseconds\n",msh.mprts.npartitions, msec/1000, msec%1000);

    for (int i = 0; i<nparts; i++){
        printf("%d %d\n",msh.mprts.parts_idx[i+1], msh.mprts.parts_idx[i+1]-msh.mprts.parts_idx[i]);
    }

    // Mesh Smoothing
    bool* bdy = Mesh_find_bdy_nodes(&msh);
    start = clock();
    Mesh_smooth2d(&msh, bdy, 500);
    diff = clock() - start;
    msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("finished mesh smoothing in %d seconds %d milliseconds\n", msec/1000, msec%1000);

    Mesh_draw(&msh);
    free(bdy);
    */

    int npoints = 2000;
    struct IntMatrix segs;
    struct DoubleMatrix xs = DoubleMatrix_create_Random(npoints,3,0.0,1.0);
    //double h = Ellipsoid(&xs, npoints,false);
    struct Mesh msh = GeoMesh_Delaunay_tet(&xs);
    bool check = check_jacobians_tet(&msh);
    Mesh2vtk_tet(&msh);
    return 0;
}