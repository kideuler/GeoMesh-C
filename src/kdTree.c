#include "GeoMesh.h"

// private functions
void kdNode_insert(struct kdNode* kdn, int nid, int node_depth, const struct DoubleMatrix* coords, int ndims);
void kdNode_free(struct kdNode* kdn);
void kdNode_find_nearest_node(struct kdNode* kdn, double* point, const struct DoubleMatrix* coords, int ndims, int* curr_nid, double* best);
void kdNode_find_target_radius(struct kdNode* kdn, double* point, const struct DoubleMatrix* coords, int ndims, int* curr_nid, double* best, double* dist, double radius, const double* h_ratios, double hgrad);

double kdTree_find_target_radius(struct kdTree* kdt, double* point){
    int ndims = kdt->coords.ncols;
    double best = 1.0e308;
    double dist = 1e308;
    int nid = -1;
    kdNode_find_target_radius(kdt->head, point, &kdt->coords, ndims, &nid, &best, &dist, kdt->radius, kdt->h_ratios, kdt->hgrad);
    return best;
}

void kdNode_find_target_radius(struct kdNode* kdn, double* point, const struct DoubleMatrix* coords, int ndims, int* curr_nid, double* best, double* dist, double radius, const double* h_ratios, double hgrad){
    if (kdn->vid < 0){ // reached end of the tree
        return;
    }
    double d, xi, r;
    if (kdn->depth == 0){
        *curr_nid = kdn->vid;
        double d = 0.0;
        for (int dim = 0; dim<ndims; dim++){
            d += pow(coords->data[ndims*kdn->vid + dim] - point[dim],2.0);
        }
        xi = d / (hgrad*radius);
        *dist = d;
        *best =  (1-(xi<1.0?xi:1.0))*h_ratios[*curr_nid]*radius + radius*(xi<1.0?xi:1.0);
    } else {
        double d = 0.0;
        for (int dim = 0; dim<ndims; dim++){
            d += pow(coords->data[ndims*kdn->vid + dim] - point[dim],2.0);
        }
        xi = d / (hgrad*radius);
        r = (1-(xi<1.0?xi:1.0))*h_ratios[*curr_nid]*radius + radius*(xi<1.0?xi:1.0);
        if (r < *best){
            *best = r;
            *dist = d;
            *curr_nid = kdn->vid;
        }
    }

    int axis = kdn->depth % ndims;
    if (point[axis] < coords->data[ndims*kdn->vid + axis]){
        kdNode_find_target_radius(kdn->left, point, coords, ndims, curr_nid, best, dist, radius, h_ratios, hgrad);
        if ((point[axis] + sqrt(*dist)) >= coords->data[ndims*kdn->vid + axis]){
            kdNode_find_target_radius(kdn->right, point, coords, ndims, curr_nid, best, dist, radius, h_ratios, hgrad);
        }
    } else {
        kdNode_find_target_radius(kdn->right, point, coords, ndims, curr_nid, best, dist, radius, h_ratios, hgrad);
        if ((point[axis] - sqrt(*dist)) <= coords->data[ndims*kdn->vid + axis]){
            kdNode_find_target_radius(kdn->left, point, coords, ndims, curr_nid, best, dist, radius, h_ratios, hgrad);
        }
    }
}

int kdTree_find_nearest_node(struct kdTree* kdt, double* point){
    int ndims = kdt->coords.ncols;
    double best = 1.0e308;
    int ans = -1;
    kdNode_find_nearest_node(kdt->head, point, &kdt->coords, ndims, &ans, &best);
    return ans;
}

void kdNode_find_nearest_node(struct kdNode* kdn, double* point, const struct DoubleMatrix* coords, int ndims, int* curr_nid, double* best){
    if (kdn->vid < 0){ // reached end of the tree
        return;
    }

    if (kdn->depth == 0){
        *curr_nid = kdn->vid;
        double d = 0.0;
        for (int dim = 0; dim<ndims; dim++){
            d += pow(coords->data[ndims*kdn->vid + dim] - point[dim],2.0);
        }
        *best = d;
    } else {
        double d = 0.0;
        for (int dim = 0; dim<ndims; dim++){
            d += pow(coords->data[ndims*kdn->vid + dim] - point[dim],2.0);
        }
        if (d < *best){
            *best = d;
            *curr_nid = kdn->vid;
        }
    }

    int axis = kdn->depth % ndims;
    if (point[axis] < coords->data[ndims*kdn->vid + axis]){
        kdNode_find_nearest_node(kdn->left, point, coords, ndims, curr_nid, best);
        if ((point[axis] + sqrt(*best)) >= coords->data[ndims*kdn->vid + axis]){
            kdNode_find_nearest_node(kdn->right, point, coords, ndims, curr_nid, best);
        }
    } else {
        kdNode_find_nearest_node(kdn->right, point, coords, ndims, curr_nid, best);
        if ((point[axis] - sqrt(*best)) <= coords->data[ndims*kdn->vid + axis]){
            kdNode_find_nearest_node(kdn->left, point, coords, ndims, curr_nid, best);
        }
    }
}


struct kdTree kdTree_create(const struct DoubleMatrix coords){
    struct kdTree kdt;
    kdt.head = (struct kdNode*) malloc(sizeof(struct kdNode));
    kdt.head->vid = -1;
    kdt.coords = coords;
    int ndims = coords.ncols;
    int nnodes = coords.nrows;
    for (int n = 0; n<nnodes; n++){
        kdNode_insert(kdt.head, n, 0, &kdt.coords, ndims);
    }

    return kdt;
}

void kdNode_insert(struct kdNode* kdn, int nid, int node_depth, const struct DoubleMatrix* coords, int ndims){
    if (kdn->vid == -1){
        kdn->vid = nid;
        kdn->depth = node_depth;
        kdn->left = (struct kdNode*) malloc(sizeof(struct kdNode));
        kdn->right = (struct kdNode*) malloc(sizeof(struct kdNode));
        kdn->left->vid = -1;
        kdn->right->vid = -1;
        //printf("depth: %d node: %d\n",node_depth,nid);
        return;
    } else {
        int axis = kdn->depth % ndims;
        if (coords->data[ndims*nid + axis] <= coords->data[ndims*kdn->vid + axis]){
            kdNode_insert(kdn->left, nid, node_depth+1, coords, ndims);
        } else {
            kdNode_insert(kdn->right, nid, node_depth+1, coords, ndims);
            return;
        }
    }
}

void kdTree_free(struct kdTree* kdt){
    kdNode_free(kdt->head);
    free(kdt->head);
}

void kdNode_free(struct kdNode* kdn){
    if (kdn->left->vid == -1){
        free(kdn->left);
    } else {
        kdNode_free(kdn->left);
        free(kdn->left);
    }

    if (kdn->right->vid == -1){
        free(kdn->right);
    } else {
        kdNode_free(kdn->right);
        free(kdn->right);
    }
    return;
}