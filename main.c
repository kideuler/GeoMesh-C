#include "GeoMesh.h"

int main(void){
    srand ( time(NULL) );
    int npoints = 40;

    struct DoubleMatrix xs;
    struct IntMatrix segs;;
    double h = Ellipse(&segs, &xs, npoints,true);
    //load_lakeSuperior(&segs, &xs);
    
    clock_t start = clock(), diff;
    //struct Mesh msh = GeoMesh_Delaunay(&xs);
    struct Mesh msh = GeoMesh_ConstrainedDelaunay(&segs,&xs);
    diff = clock() - start;
    int msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("finished initial triangulation\n");
    printf("Created %d Triangles in %d seconds %d milliseconds\n", msh.nelems, msec/1000, msec%1000);


    bool is_sib = check_sibhfs(&msh);

    start = clock();
    GeoMesh_DelaunayRefine(&msh,true, h,1);
    diff = clock() - start;
    msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("finished Delaunay refinement\n");
    printf("Created %d Triangles in %d seconds %d milliseconds\n", msh.nelems, msec/1000, msec%1000);

    Mesh_draw(&msh);
    return 0;
}