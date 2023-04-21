#include "GeoMesh.h"

struct DoubleMatrix Ellipse(int npoints){
    struct DoubleMatrix xs = DoubleMatrix_create(npoints,2);
    double t;
    for (int i = 0; i<npoints; i++){
        t = 2*M_PI*(((double) i) / ((double) npoints));
        xs.data[2*i] = 0.2*cos(t)+0.5;
        xs.data[2*i+1] = 0.1*sin(t)+0.5;
    }

    return xs;
}

struct DoubleMatrix Flower(int npoints){
    struct DoubleMatrix xs = DoubleMatrix_create(npoints,2);
    double t;
    for (int i = 0; i<npoints; i++){
        t = 2*M_PI*(((double) i) / ((double) npoints));
        xs.data[2*i] = (0.25 + 0.1*sin(5*t))*cos(t)/2.0+0.5;
        xs.data[2*i+1] = (0.25 + 0.1*sin(5*t))*sin(t)/2.0+0.5;
    }

    return xs;
}

int main(void){
    srand ( time(NULL) );
    int npoints = 100;
    //struct DoubleMatrix xs = DoubleMatrix_create_Random(npoints,2,0.3,0.7);
    struct DoubleMatrix xs = Flower(npoints);
    struct IntMatrix segs = IntMatrix_create(npoints,2);
    for(int i = 0; i<npoints; i++){segs.data[2*i]=i; segs.data[2*i+1]=(i+1)%npoints;}
    clock_t start = clock(), diff;
    //struct Mesh msh = GeoMesh_Delaunay(&xs);
    struct Mesh msh = GeoMesh_ConstrainedDelaunay(&segs,&xs);
    diff = clock() - start;
    int msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("finished initial triangulation\n");
    printf("Created %d Triangles in %d seconds %d milliseconds\n", msh.nelems, msec/1000, msec%1000);
    Mesh_compute_OneringElements(&msh, 15);
    //IntMatrix_print(&msh.stncl.hvids);
    Mesh_draw(&msh);
    return 0;
}