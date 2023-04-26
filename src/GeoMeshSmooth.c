#include "GeoMesh.h"

void Mesh_smooth2d(struct Mesh* msh, bool* no_move, int max_Niters){
    if (!msh->hasStencil){
        Mesh_compute_OneringElements(msh, 10);
    }

    int iter = 0;
    double Energy = 1e12;
    double Energy_old = 1e14;
    double det,alpha;
    int iter_move;

    int nv = msh->coords.nrows;

    // create grads arrays on the heap
    double* Grads = malloc(2*nv*sizeof(double));
    memset(Grads, 0.0, 2*nv*sizeof(double));
    double* Grad_elem = malloc(3*2*sizeof(double));

    // create hessian arrays on the heap
    double* Hessian = malloc(2*2*nv*sizeof(double));
    memset(Hessian, 0.0, 2*2*nv*sizeof(double));
    double* Hess_elem = malloc(3*2*2*sizeof(double));
    double* Hess_inv = malloc(2*2*sizeof(double));

    // create coords_diff on the heap
    double* coords_diff = malloc(nv*2*sizeof(double));

    double ps[6];

    while (iter < max_Niters && Energy_old > Energy){
        // iteration of Energy based smoothing
        double Energy_total = 0.0;
        double Energy;
        int v;

        // Accumulating Energy
        for (int n = 0; n<msh->nelems; n++){
            ps[0] = msh->coords.data[msh->elems.data[3*n]];
            ps[1] = msh->coords.data[msh->elems.data[3*n]+1];
            ps[2] = msh->coords.data[msh->elems.data[3*n+1]];
            ps[3] = msh->coords.data[msh->elems.data[3*n+1]+1];
            ps[4] = msh->coords.data[msh->elems.data[3*n+2]];
            ps[5] = msh->coords.data[msh->elems.data[3*n+2]+1];

            Energy = isometry_energy_tri(ps, Grad_elem, Hess_elem);
            for (int ii = 0; ii<3; ii++){
                v = msh->elems.data[3*n+ii];
                Grads[v] += Grad_elem[ii];
                Grads[nv + v] += Grad_elem[3+ii];

                Hessian[4*nv] += Hess_elem[4*ii];
                Hessian[4*nv+1] += Hess_elem[4*ii+1];
                Hessian[4*nv+2] += Hess_elem[4*ii+2];
                Hessian[4*nv+3] += Hess_elem[4*ii+3];
            }
            Energy_total += Energy;
        }

        // Moving nodes
        for (int n = 0; n<nv; n++){
            if (!no_move[n]){
                det = Hessian[4*n]*Hessian[4*n+3] - Hessian[4*n+1]*Hessian[4*n+2];
                if (abs(det) > 1e-3){
                    Hess_inv[0] = Hessian[4*n+3]/det;
                    Hess_inv[1] = -Hessian[4*n+1]/det;
                    Hess_inv[2] = -Hessian[4*n+2]/det;
                    Hess_inv[3] = Hessian[4*n]/det;

                    coords_diff[2*n] = -Hess_inv[0]*Grads[n] - Hess_inv[1]*Grads[n + nv];
                }
            }
        }


        memset(Grads, 0.0, 2*nv*sizeof(double));
        memset(Hessian, 0.0, 2*2*nv*sizeof(double));
    }


    free(coords_diff);
    free(Grads); free(Grad_elem);
    free(Hessian); free(Hess_elem);
    free(Hess_inv);
}

double isometry_energy_tri(const double ps[6], double* Grad_elem, double* Hess_elem){
    double Energy;

    return Energy;
}