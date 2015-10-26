#include "header.h"

Parameters *set_default_params(void)
{
  int i;
  double *arr;
  Parameters *params = malloc(sizeof(Parameters));

  params->n_grid = 20;
  params->h = 0.5;
  params->L = params->h*params->n_grid;
  parmas->G = 2;
  params->macro_xs_f = calloc(params->G, sizeof(double));
  params->macro_xs_a = calloc(params->G, sizeof(double));
  params->macro_xs_e = malloc(params->G*sizeof(double*));
  arr = calloc(params->G*params->G, sizeof(double));
  for(i=0; i<params->G; i++){
    params->macro_xs_e[i] = arr;
    arr += params->G;
  }
  params->macro_xs_t = calloc(params->G, sizeof(double));
  params->macro_xs_r = calloc(params->G, sizeof(double));
  params->mu = 0.0;
  params->nu = calloc(params->G, sizeof(double));
  params->k = 1.0;
  params->D = calloc(params->G, sizeof(double));
  params->chi = calloc(params->G, sizeof(double));
  params->bc = VACUUM;
  params->max_inner = 1000;
  params->max_outer = 1000;
  params->thresh = 0.00001;
  params->write_flux = 0;
  params->flux_file = NULL;

  return params;
}

// Initial guess of uniform flux
double ****init_flux(Parameters *params)
{
  int g, i, j, k;
  double ****phi;

  phi = matrix4D(params->G, params->n_grid, params->n_grid, params->n_grid);

  for(g=0; g<params->G; g++){
    for(i=0; i<params->n_grid; i++){
      for(j=0; j<params->n_grid; j++){
        for(k=0; k<params->n_grid; k++){
          phi[g][i][j][k] = 1.0;
        }
      }
    }
  }

  return phi;
}

void free_params(Parameters *params)
{
  free(params->macro_xs_a);
  free(params->macro_xs_e[0]);
  free(params->macro_xs_e);
  free(params->macro_xs_f);
  free(params->macro_xs_t);
  free(params->macro_xs_r);
  free(params->D);
  free(params->nu);
  free(params->chi);
  free(params);
}

void free_flux(double ***phi)
{
  free_matrix4D(phi);
}
