#include "header.h"

Parameters *set_default_params(void)
{
  Parameters *params = malloc(sizeof(Parameters));

  params->n_grid = 20;
  params->h = 0.5;
  params->L = params->h*params->n_grid;
  params->xs_f = 2.29;
  params->xs_a = 3.42;
  params->xs_s = 2.29;
  params->xs_t = params->xs_f + params->xs_a + params->xs_s;
  params->mu = 0.0;
  params->nu = 1.5;
  params->k = 1.0;
  params->D = 1/(3*params->xs_t - params->mu*params->xs_s);
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
