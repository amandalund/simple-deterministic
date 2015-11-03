#include "header.h"

Parameters *set_default_params(void)
{
  Parameters *params = malloc(sizeof(Parameters));

  params->n_grid = 20;
  params->h = 0.5;
  params->L = params->h*params->n_grid;
  params->G = 1;
  params->mu = 0.0;
  params->k = 1.0;
  params->bc = VACUUM;
  params->max_inner = 1000;
  params->max_outer = 1000;
  params->thresh = 0.00001;
  params->write_flux = 0;
  params->flux_file = NULL;
  params->group_file = NULL;

  return params;
}

// Initial guess of uniform flux
double ****init_flux(Parameters *params)
{
  int g, i, j, k;
  double ****phi;
  double val = 1.0/params->G;

  phi = matrix4D(params->G, params->n_grid, params->n_grid, params->n_grid);

  for(g=0; g<params->G; g++){
    for(i=0; i<params->n_grid; i++){
      for(j=0; j<params->n_grid; j++){
        for(k=0; k<params->n_grid; k++){
          phi[g][i][j][k] = val;
        }
      }
    }
  }

  return phi;
}

void free_params(Parameters *params)
{
  free(params->xs_a);
  free(params->xs_s[0]);
  free(params->xs_s);
  free(params->xs_f);
  free(params->xs_t);
  free(params->xs_r);
  free(params->D);
  free(params->nu);
  free(params->chi);
  free(params);
}

void free_flux(double ****phi)
{
  free_matrix4D(phi);
}
