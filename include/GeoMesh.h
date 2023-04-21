#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <pbPlots.h>
#include <supportLib.h>

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

struct Stencil {
    int maxne;
    struct IntMatrix hvids; // encodes element and local vertex id
    int* nhvids; // number of elements for each node
};

struct Mesh {
    int nelems;
    struct DoubleMatrix coords;
    struct IntMatrix elems;
    struct IntMatrix sibhfs; // AHF data structure
    bool *delete_elem;
    bool *bwork;
    struct Stencil stncl;
};

// matrix functions
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
struct Mesh GeoMesh_ConstrainedDelaunay(struct IntMatrix *segments, struct DoubleMatrix *xs);

// Delaunay subfunctions
void Mesh_compute_OneringElements(struct Mesh* msh, int maxne);
void Mesh_deleteElems(struct Mesh* msh);
void Recursive_findDelete(struct Mesh* msh, int hfid);
void Mesh_flip_insertion(struct Mesh* msh, int* vid, int tri_start);
bool Mesh_find_enclosing_tri(struct Mesh* msh, int* tri, double ps[2]);
void flip_edge(struct Mesh* msh, int eid, int lid);
bool inside_tri(const double xs[3][2], const double ps[2]);
bool inside_circumtri(const double xs[3][2], const double ps[2]);
double* circumcenter(const double xs[3][2]);
static bool Line_cross(double p1x, double p1y, double p2x, double p2y, \
double p3x, double p3y, double p4x, double p4y);

// drawing mesh png using pbPlot
void Mesh_draw(struct Mesh* msh);

#endif