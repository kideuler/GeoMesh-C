#include "GeoMesh.h"

void Mesh_smooth2d(struct Mesh* msh, bool* no_move, int max_Niters){
    if (!msh->hasStencil){
        Mesh_compute_OneringElements(msh, 10);
    }

    int iter = 0;
    int v;
    double Energy, Energy_total;
    double Energy_old = 1e14;
    double det,alpha;
    int iter_move;

    int nv = msh->coords.nrows;

    // create grads arrays on the heap
    double* Grads = (double*) malloc(2*nv*sizeof(double));
    memset(Grads, 0.0, 2*nv*sizeof(double));
    double* Grad_elem = (double*) malloc(3*2*sizeof(double));
    memset(Grad_elem, 0.0, 3*2*sizeof(double));

    // create hessian arrays on the heap
    double* Hessian = (double*) malloc(2*2*nv*sizeof(double));
    memset(Hessian, 0.0, 2*2*nv*sizeof(double));
    double* Hess_elem = (double*) malloc(3*2*2*sizeof(double));
    memset(Hess_elem, 0.0, 2*2*3*sizeof(double));
    double* Hess_inv = (double*) malloc(2*2*sizeof(double));

    // create coords_diff on the heap
    double* coords_diff = (double*) malloc(nv*2*sizeof(double));
    memset(coords_diff, 0.0, nv*2*sizeof(double));

    double ps[6];

    while (iter < max_Niters){
        // iteration of Energy based smoothing
        Energy_total = 0.0;
        
        // Accumulating Energy
        for (int n = 0; n<msh->nelems; n++){
            ps[0] = msh->coords.data[2*msh->elems.data[3*n]];
            ps[1] = msh->coords.data[2*msh->elems.data[3*n]+1];
            ps[2] = msh->coords.data[2*msh->elems.data[3*n+1]];
            ps[3] = msh->coords.data[2*msh->elems.data[3*n+1]+1];
            ps[4] = msh->coords.data[2*msh->elems.data[3*n+2]];
            ps[5] = msh->coords.data[2*msh->elems.data[3*n+2]+1];

            Energy = isometry_energy_tri(ps, Grad_elem, Hess_elem);
            for (int ii = 0; ii<3; ii++){
                v = msh->elems.data[3*n+ii];
                Grads[v] += Grad_elem[ii];
                Grads[nv + v] += Grad_elem[3+ii];

                Hessian[4*v] += Hess_elem[4*ii];
                Hessian[4*v+1] += Hess_elem[4*ii+1];
                Hessian[4*v+2] += Hess_elem[4*ii+2];
                Hessian[4*v+3] += Hess_elem[4*ii+3];
            }
            Energy_total += Energy;
        }

        // check if loop should end
        if (Energy_total >= Energy_old){
            if (Energy_total > Energy_old){
                for (int n = 0; n<nv; n++){
                    msh->coords.data[2*n] -= coords_diff[2*n];
                    msh->coords.data[2*n+1] -= coords_diff[2*n+1];
                }
            }
            break;
        } else {
            //printf("Isometry Energy: %f\n", Energy_total);
        }
        Energy_old = Energy_total;

        // Calcing differentials
        for (int n = 0; n<nv; n++){
            if (!no_move[n]){
                det = Hessian[4*n]*Hessian[4*n+3] - Hessian[4*n+1]*Hessian[4*n+2];
                det = (det<0)?(-det):(det);
                if (det > 1e-3){
                    Hess_inv[0] = Hessian[4*n+3]/det;
                    Hess_inv[1] = -Hessian[4*n+1]/det;
                    Hess_inv[2] = -Hessian[4*n+2]/det;
                    Hess_inv[3] = Hessian[4*n]/det;

                    coords_diff[2*n] = -Hess_inv[0]*Grads[n] - Hess_inv[1]*Grads[n + nv];
                    coords_diff[2*n+1] = -Hess_inv[2]*Grads[n] - Hess_inv[3]*Grads[n + nv];
                    alpha = 0.1;
                    iter_move = 0;
                    while (false){
                        alpha = 0.5*alpha;
                        iter_move++;
                    }

                    if (iter_move < 10){
                        coords_diff[2*n] = coords_diff[2*n]*alpha;
                        coords_diff[2*n+1] = coords_diff[2*n+1]*alpha;
                    }
                }
            }
        }

        // moving nodes
        for (int n = 0; n<nv; n++){
            msh->coords.data[2*n] += coords_diff[2*n];
            msh->coords.data[2*n+1] += coords_diff[2*n+1];
        }

        memset(Grads, 0.0, 2*nv*sizeof(double));
        memset(Hessian, 0.0, 2*2*nv*sizeof(double));
        iter++;
    }


    free(coords_diff);
    free(Grads); free(Grad_elem);
    free(Hessian); free(Hess_elem);
    free(Hess_inv);
}

double isometry_energy_tri(const double xs[6], double* Grad, double* Hess){
    double Energy;

    double e12[2],e12_orth[2],e23[2],e23_orth[2],e31[2],e31_orth[2], sql12, sql23, sql31;
    e12[0] = xs[2]-xs[0]; e12[1] = xs[3]-xs[1];
    sql12 = e12[0]*e12[0] + e12[1]*e12[1];
    e23[0] = xs[4]-xs[2]; e23[1] = xs[5]-xs[3];
    sql23 = e23[0]*e23[0] + e23[1]*e23[1];
    e31[0] = xs[0]-xs[4]; e31[1] = xs[1]-xs[5];
    sql31 = e31[0]*e31[0] + e31[1]*e31[1];
    double area2 = (e12[0]*e23[1] - e12[1]*e23[0]);
    double area = 0.5*area2;
    e12_orth[0] = -e12[1]; e12_orth[1] = e12[0];
    e23_orth[0] = -e23[1]; e23_orth[1] = e23[0];
    e31_orth[0] = -e31[1]; e31_orth[1] = e31[0];
    double cts[3] = {1/(sqrt(3.0)*area),1/(sqrt(3.0)*area),1/(sqrt(3.0)*area)};

    Energy = (cts[2]*sql12 + cts[0]*sql23 + cts[1]*sql31);

    int i,j;
    for (i=0; i<2; i++){
        for(j=0; j<3; j++){
            Grad[3*i + j] = 0.0;
        } 
    }

    Grad[0] = -2*cts[2]*e12[0];
    Grad[3] = -2*cts[2]*e12[1];
    Grad[1] = 2*cts[2]*e12[0];
    Grad[4] = 2*cts[2]*e12[0];

    Grad[1] += -2*cts[0]*e23[0];
    Grad[4] += -2*cts[0]*e23[1];
    Grad[2] = 2*cts[0]*e23[0];
    Grad[5] = 2*cts[0]*e23[1];

    Grad[2] += -2*cts[1]*e31[0]; 
    Grad[5] += -2*cts[1]*e31[1];
    Grad[0] += 2*cts[1]*e31[0];
    Grad[3] += 2*cts[1]*e31[1];
    double energy_a = Energy/area2;

    Grad[0] -= energy_a*e23_orth[0];
    Grad[3] -= energy_a*e23_orth[1];
    Grad[1] -= energy_a*e31_orth[0];
    Grad[4] -= energy_a*e31_orth[1];
    Grad[2] -= energy_a*e12_orth[0];
    Grad[5] -= energy_a*e12_orth[1];

    double c = 2*(cts[2]+cts[1]);
    double C[4];
    C[0] = c; C[1] = 0.0; C[2] = 0.0; C[3] = c;

    for (i=0;i<2;i++){
        for(j=0;j<2;j++){
            Hess[2*i + j] = C[2*i + j] - (Grad[3*i]*e23_orth[j]/area2 + Grad[3*j]*e23_orth[i]/area2);
        }
    }
    for (i = 0; i<2; i++){
        for (j = 0; j<2; j++){
            Hess[2*i + j + 4] = C[2*i + j] - (Grad[3*i+1]*e31_orth[j]/area2 + Grad[3*j+1]*e31_orth[i]/area2);
        }
    }
    for (i = 0; i<2; i++){
        for (j = 0; j<2; j++){
            Hess[2*i + j + 8] = C[2*i+j] - (Grad[3*i+2]*e12_orth[j]/area2 + Grad[3*j+2]*e12_orth[i]/area2);
        }
    }

    return Energy;
}