#include "header.h"

// Gauss-Seidel iterative solver
void solve(double ****phi, Parameters *params)
{
  double ***S, ***S0;             // current and previous estimates of source
  double S_norm, S0_norm;         // current and previous source vector norm (||phi||_1)
  double S_max;                   // source normalization (||S||_inf = max(S))
  double k_prev;                  // previous estimate of keff
  double err_outer;               // error to determine when to break outer iterations
  int i, j, k;                    // indices over flux elements
  int l;                          // index over outer iterations
  int g;                          // index over energy group
  int n;                          // total number of mesh elements

  n = params->n_grid*params->n_grid*params->n_grid;
  S = matrix3D(params->n_grid, params->n_grid, params->n_grid);
  S0 = matrix3D(params->n_grid, params->n_grid, params->n_grid);
  compute_source(phi, S0, params, &S_max);

  for(l=0; l<params->max_outer; l++){

    err_outer = S_norm = S0_norm = 0;

    // Solve M*phi^(n+1) = (1/k)*S^n for each group
    for(g=0; g<params->G; g++){
      solve_inner(phi, S0, params, g);
    }

    // Calculate the new fission source
    compute_source(phi, S, params, &S_max);

    // Calculate source norms, error, and normalize source
    for(i=0; i<params->n_grid; i++){
      for(j=0; j<params->n_grid; j++){
        for(k=0; k<params->n_grid; k++){

          S_norm += S[i][j][k];
          S0_norm += S0[i][j][k];
          err_outer += fabs(S[i][j][k] - S0[i][j][k]);

          // Normalize source
          S[i][j][k] /= S_max;
          for(g=0; g<params->G; g++){
            phi[g][i][j][k] /= S_max;
          }
        }
      }
    }

    k_prev = params->k;
    params->k = k_prev*S_norm/S0_norm;
    memcpy(S0[0][0], S[0][0], n*sizeof(double));

    if(err_outer/n < params->thresh && fabs(params->k - k_prev)/params->k < params->thresh) break;
    if(l == params->max_outer - 1) printf("WARNING: did not converge in %d outer iterations\n", params->max_outer);
  }

  free_matrix3D(S);
  free_matrix3D(S0);
}

void solve_inner(double ****phi, double ***S, Parameters *params, int g)
{
  double c, d;                    // matrix coefficients
  double phi_prev;                // flux element from previous iteration
  double err_inner;               // error to determine when to break inner iterations
  double t1, t2, t3, t4, t5, t6;  // neighboring mesh elements
  double rhs;                     // rhs term
  int i, j, k;                    // indices over flux elements
  int f;                          // index over previous groups
  int m;                          // index over inner iterations
  int n;                          // total number of mesh elements

  n = params->n_grid*params->n_grid*params->n_grid;
  c = params->D[g]/(params->h*params->h);
  d = 1/(6*c + params->macro_xs_r[g]);
//  d = 1/(6*c + params->macro_xs_a[g]);

  for(m=0; m<params->max_inner; m++){

    err_inner = 0;

    for(i=0; i<params->n_grid; i++){
      for(j=0; j<params->n_grid; j++){
        for(k=0; k<params->n_grid; k++){

          phi_prev = phi[g][i][j][k];

          t1 = t2 = t3 = t4 = t5 = t6 = FLAG;

          // boundary conditions
          if(i == params->n_grid-1){
            if(params->bc == VACUUM) t1 = 0;
            else if(params->bc == PERIODIC) t1 = phi[g][0][j][k];
          }
          if(i == 0){
            if(params->bc == VACUUM) t2 = 0;
            else if(params->bc == PERIODIC) t2 = phi[g][params->n_grid-1][j][k];
          }
          if(j == params->n_grid-1){
            if(params->bc == VACUUM) t3 = 0;
            else if(params->bc == PERIODIC) t3 = phi[g][i][0][k];
          }
          if(j == 0){
            if(params->bc == VACUUM) t4 = 0;
            else if(params->bc == PERIODIC) t4 = phi[g][i][params->n_grid-1][k];
          }
          if(k == params->n_grid-1){
            if(params->bc == VACUUM) t5 = 0;
            else if(params->bc == PERIODIC) t5 = phi[g][i][j][0];
          }
          if(k == 0){
            if(params->bc == VACUUM) t6 = 0;
            else if(params->bc == PERIODIC) t6 = phi[g][i][j][params->n_grid-1];
          }

          if(t1 == FLAG) t1 = phi[g][i+1][j][k];
          if(t2 == FLAG) t2 = phi[g][i-1][j][k];
          if(t3 == FLAG) t3 = phi[g][i][j+1][k];
          if(t4 == FLAG) t4 = phi[g][i][j-1][k];
          if(t5 == FLAG) t5 = phi[g][i][j][k+1];
          if(t6 == FLAG) t6 = phi[g][i][j][k-1];

          rhs = params->nu[g]*params->macro_xs_f[g]*params->chi[g]/params->k*S[i][j][k];

          // Add contribution from neutrons scattering in from higher groups
          for(f=0; f<g; f++){
            rhs += params->macro_xs_e[f][g]*phi[f][i][j][k];
          }

          phi[g][i][j][k] = c*d*(t1 + t2 + t3 + t4 + t5 + t6) + d*rhs;

          err_inner += fabs(phi[g][i][j][k] - phi_prev);
        }
      }
    }

    if(err_inner/n < params->thresh) break;
    if(m == params->max_inner - 1) printf("WARNING: did not converge in %d inner iterations\n", params->max_inner);
  }

  return;
}

void compute_source(double ****phi, double ***S, Parameters *params, double *max)
{
  int g, i, j, k;

  *max = 0;

  // Accumulate source from each group's flux
  for(i=0; i<params->n_grid; i++){
    for(j=0; j<params->n_grid; j++){
      for(k=0; k<params->n_grid; k++){
        S[i][j][k] = 0;
        for(g=0; g<params->G; g++){
          S[i][j][k] += phi[g][i][j][k];
        }
        if(S[i][j][k] > *max) *max = S[i][j][k];
      }
    }
  }

  return;
}
