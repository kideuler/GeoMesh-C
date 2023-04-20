#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#ifndef GEOMESH
#define GEOMESH

struct DoubleMatrix{
    int nrows;
    int ncols;
    double* data;
};

struct IntMatrix{
    int nrows;
    int ncols;
    int* data;
};

struct Mesh {
    int nelems;
    struct DoubleMatrix coords;
    struct IntMatrix elems;
    struct IntMatrix sibhfs;
    bool *delete_elem;
};


struct DoubleMatrix DoubleMatrix_create(int nrows, int ncols);
struct IntMatrix IntMatrix_create(int nrows, int ncols);

double drand(double low, double high);
struct DoubleMatrix DoubleMatrix_create_Random(int nrows, int ncols, double low, double high);

void DoubleMatrix_print(struct DoubleMatrix* M);
void IntMatrix_print(struct IntMatrix* M);

double DoubleMatrix_min(struct DoubleMatrix* M, int col);
double DoubleMatrix_max(struct DoubleMatrix* M, int col);

// printing mesh
void Mesh_print(struct Mesh *msh);

// ahf functions
int elids2hfid(int eid, int lid);
int hfid2eid(int hfid);
int hfid2lid(int hfid);

// Delaunay triangulation
struct Mesh GeoMesh_Delaunay(struct DoubleMatrix *xs);
void Mesh_flip_insertion(struct Mesh* msh, int* vid, int tri_start);
bool Mesh_find_enclosing_tri(struct Mesh* msh, int* tri, double ps[2]);
bool inside_tri(const double xs[3][2], const double ps[2]);

#endif