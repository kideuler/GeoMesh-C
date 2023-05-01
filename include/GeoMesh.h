#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <pbPlots.h>
#include <supportLib.h>
#include <assert.h>
//#include <MagickCore/MagickCore.h>

#ifndef GEOMESH
#define GEOMESH

// structure for double matrix
struct DoubleMatrix{
    int nrows;
    int ncols;
    double* data;
};

// structure for matrix of ints
struct IntMatrix{
    int nrows;
    int ncols;
    int* data;
};

struct kdNode {
    int depth; // determines the axis used
    int vid; // node id of this kdnode
    struct kdNode *left; // left pointer
    struct kdNode *right; // right pointer
};

struct kdTree { // data structure designed for quick target circumcircle target radius

    struct kdNode* head; // head node of the tree
    struct DoubleMatrix coords; // coordinate array
    double h; // average h of the mesh
    double* h_ratios; // h fraction of each coordinate
    double hgrad; // gradiation coefficient 
    
    // HOW TO CALCULATE radius_target:
    // r = distance from coords[vid] to query node
    // xi = r / (hgrad*h) local variable xi
    // radius_target = 1-min(xi,1.0)*h_ratios[vid]*h + h*min(xi,1.0);
};

// stencil used to ensure 
struct Stencil {
    int maxne;
    struct IntMatrix hvids; // encodes element and local vertex id
    int* nhvids; // number of elements for each node
};

struct Mesh {
    int nelems;
    bool hasStencil;
    struct DoubleMatrix coords;
    struct IntMatrix elems;
    struct IntMatrix sibhfs; // AHF data structure
    bool *delete_elem;
    bool *bwork;
    bool *on_boundary;
    int* stack; // stack buffer for insertion
    struct Stencil stncl;
};

// Example functions to generate points and segments
double Ellipse(struct IntMatrix* segments, struct DoubleMatrix* coords, int npoints, bool box);
double Flower(struct IntMatrix* segments, struct DoubleMatrix* coords, int npoints, bool box);
double Airfoil(struct IntMatrix* segments, struct DoubleMatrix* coords, int npoints, bool box);

// matrix functions
struct DoubleMatrix DoubleMatrix_create(int nrows, int ncols);
struct IntMatrix IntMatrix_create(int nrows, int ncols);
void DoubleMatrix_resize(struct DoubleMatrix* M, int newsz);
void IntMatrix_resize(struct IntMatrix* M, int newsz);

double drand(double low, double high);
struct DoubleMatrix DoubleMatrix_create_Random(int nrows, int ncols, double low, double high);

void DoubleMatrix_print(struct DoubleMatrix* M);
void IntMatrix_print(struct IntMatrix* M);

double DoubleMatrix_min(struct DoubleMatrix* M, int col);
double DoubleMatrix_max(struct DoubleMatrix* M, int col);

// resizing mesh
void Mesh_resize(struct Mesh* msh, int newsz);

// printing mesh
void Mesh_print(struct Mesh *msh);

// ahf functions
int elids2hfid(int eid, int lid);
int hfid2eid(int hfid);
int hfid2lid(int hfid);

// Delaunay triangulation
struct Mesh GeoMesh_Delaunay(struct DoubleMatrix *xs);
struct Mesh GeoMesh_ConstrainedDelaunay(struct IntMatrix *segments, struct DoubleMatrix *xs);
void GeoMesh_DelaunayRefine(struct Mesh *msh, bool use_edgelengh, double h_target, int point_algorithm);

// isoparametric mesh smoothing functions
void Mesh_smooth2d(struct Mesh* msh, bool* no_move, int max_Niters);
double isometry_energy_tri(const double ps[6], double* Grad_elem, double* Hess_elem);

// Delaunay subfunctions
void Mesh_compute_OneringElements(struct Mesh* msh, int maxne);
void Mesh_deleteElems(struct Mesh* msh);
void Recursive_findDelete(struct Mesh* msh, int hfid);
void Mesh_flip_insertion(struct Mesh* msh, int* vid, int tri_start);
void Mesh_flip_insertion_segment(struct Mesh* msh, int vid, int hfid);
bool Mesh_find_enclosing_tri(struct Mesh* msh, int* tri, double ps[2]);
void flip_edge(struct Mesh* msh, int eid, int lid);
bool inside_tri(const double xs[3][2], const double ps[2]);
bool inside_circumtri(const double xs[3][2], const double ps[2]);
bool inside_diametral(struct Mesh* msh, int hfid, double ps[2]);
void circumcenter(const double xs[3][2], double* C);
void off_circumcenter(const double xs[3][2], double beta, double* C);
double area_tri(double xs[3][2]);
static bool Line_cross(double p1x, double p1y, double p2x, double p2y, \
double p3x, double p3y, double p4x, double p4y);
bool convex_quad(struct Mesh* msh, int eid, int lid);
double eval_alpha(const double xs[3][2],double radius_target);
double eval_trishape(const double xs[3][2]);

// mesh utility functions
bool check_sibhfs(struct Mesh* msh);
bool check_jacobians(struct Mesh* msh);
bool* Mesh_find_bdy_nodes(struct Mesh* msh);

// drawing mesh png using pbPlot
void Mesh_draw(struct Mesh* msh);

// loading lake superior data into matrices
void load_lakeSuperior(struct IntMatrix* segments, struct DoubleMatrix* coords);
#endif