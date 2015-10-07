#include "header.h"

Parameters *set_default_params(void)
{
  Parameters *params = malloc(sizeof(Parameters));

  params->n_grid = 20;
  params->h = 0.5;
  params->L = params->h*params->n_grid;
  params->macro_xs_f = 2.29;
  params->macro_xs_a = 3.42;
  params->macro_xs_e = 2.29;
  params->macro_xs_t = params->macro_xs_f + params->macro_xs_a + params->macro_xs_e;
  params->mu = 0.0;
  params->nu = 1.5;
  params->k = 1.0;
  params->D = 1/(3*params->macro_xs_t - params->mu*params->macro_xs_e);
  params->bc = VACUUM;
  params->max_inner = 1000;
  params->max_outer = 1000;
  params->thresh = 0.00001;
  params->write_flux = 0;
  params->flux_file = NULL;

  return params;
}

// Initial guess of uniform flux
double ***init_flux(Parameters *params)
{
  int i, j, k;
  double ***phi;

  phi = matrix3D(params->n_grid, params->n_grid, params->n_grid);

  for(i=0; i<params->n_grid; i++){
    for(j=0; j<params->n_grid; j++){
      for(k=0; k<params->n_grid; k++){
        phi[i][j][k] = 1.0;
      }
    }
  }

  return phi;
}

void free_flux(double ***phi)
{
  free_matrix3D(phi);
}
