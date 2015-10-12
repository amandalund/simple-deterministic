#include "header.h"

// Gauss-Seidel iterative solver
void solve(double ***phi0, Parameters *params)
{
  double c, d;
  double ***phi;
  double phi_prev, k_prev;
  double err_inner, err_outer;
  double s, s0;                   // current and previous flux vector norm (||phi||_1)
  double t1, t2, t3, t4, t5, t6;  // neighboring mesh elements
  double norm;                    // flux normalization (||phi||_inf = max(phi))
  int i, j, k;                    // indices over flux elements
  int l, m;                       // indices over outer iterations, inner iterations
  int n;                          // total number of mesh elements

  c = params->D/(params->h*params->h);
  d = 1/(6*c + params->macro_xs_a);
  n = params->n_grid*params->n_grid*params->n_grid;
  phi = matrix3D(params->n_grid, params->n_grid, params->n_grid);

  for(l=0; l<params->max_outer; l++){

    err_outer = s = s0 = 0;
    memcpy(phi[0][0], phi0[0][0], n*sizeof(double));

    for(m=0; m<params->max_inner; m++){

      err_inner = norm = 0;

      for(i=0; i<params->n_grid; i++){
        for(j=0; j<params->n_grid; j++){
          for(k=0; k<params->n_grid; k++){

            phi_prev = phi[i][j][k];

            t1 = t2 = t3 = t4 = t5 = t6 = FLAG;

            // boundary conditions
            if(i == params->n_grid-1){
              if(params->bc == VACUUM) t1 = 0;
              else if(params->bc == PERIODIC) t1 = phi[0][j][k];
            }
            if(i == 0){
              if(params->bc == VACUUM) t2 = 0;
              else if(params->bc == PERIODIC) t2 = phi[params->n_grid-1][j][k];
            }
            if(j == params->n_grid-1){
              if(params->bc == VACUUM) t3 = 0;
              else if(params->bc == PERIODIC) t3 = phi[i][0][k];
            }
            if(j == 0){
              if(params->bc == VACUUM) t4 = 0;
              else if(params->bc == PERIODIC) t4 = phi[i][params->n_grid-1][k];
            }
            if(k == params->n_grid-1){
              if(params->bc == VACUUM) t5 = 0;
              else if(params->bc == PERIODIC) t5 = phi[i][j][0];
            }
            if(k == 0){
              if(params->bc == VACUUM) t6 = 0;
              else if(params->bc == PERIODIC) t6 = phi[i][j][params->n_grid-1];
            }

            if(t1 == FLAG) t1 = phi[i+1][j][k];
            if(t2 == FLAG) t2 = phi[i-1][j][k];
            if(t3 == FLAG) t3 = phi[i][j+1][k];
            if(t4 == FLAG) t4 = phi[i][j-1][k];
            if(t5 == FLAG) t5 = phi[i][j][k+1];
            if(t6 == FLAG) t6 = phi[i][j][k-1];

            phi[i][j][k] = c*d*(t1 + t2 + t3 + t4 + t5 + t6) + d*params->nu*params->macro_xs_f/params->k*phi0[i][j][k];

            err_inner += fabs(phi[i][j][k] - phi_prev);

            if(phi[i][j][k] > norm) norm = phi[i][j][k];
          }
        }
      }

      if(err_inner/n < params->thresh) break;
      if(m == params->max_inner - 1) printf("WARNING: did not converge in %d inner iterations\n", params->max_inner);

    }

    for(i=0; i<params->n_grid; i++){
      for(j=0; j<params->n_grid; j++){
        for(k=0; k<params->n_grid; k++){

          s += phi[i][j][k];
          s0 += phi0[i][j][k];
          err_outer += fabs(phi[i][j][k] - phi0[i][j][k]);

          // Normalize flux
          phi[i][j][k] /= norm;
        }
      }
    }

    k_prev = params->k;
    params->k = k_prev*s/s0;
    memcpy(phi0[0][0], phi[0][0], n*sizeof(double));
    if(err_outer/n < params->thresh && fabs(params->k - k_prev)/params->k < params->thresh) break;
    if(l == params->max_outer - 1) printf("WARNING: did not converge in %d outer iterations\n", params->max_outer);
  }

  free_matrix3D(phi);
}
