#include "GeoMesh.h"

int main(void){
    struct DoubleMatrix xs = DoubleMatrix_create_Random(10,2,-10.0,10.0);
    struct Mesh msh = GeoMesh_Delaunay(&xs);
    Mesh_print(&msh);
}