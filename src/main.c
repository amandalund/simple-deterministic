#include "header.h"

int main(int argc, char *argv[])
{
  Parameters *params; // user defined parameters
  double ****phi;     // flux array
  FILE *fp = NULL;    // output file

  // Get inputs
  params = set_default_params();
  parse_params("parameters", params);
  read_CLI(argc, argv, params);
  print_params(params);

  // Initial guess of flux
  phi = init_flux(params);

  solve(phi, params);

  printf("keff = %f\n", params->k);

  // Write solution
  if(params->write_flux == TRUE){
    write_flux(phi, params, fp);
  }

  // Free memory
  free_flux(phi);
  free_params(params);

  return 0;
}
