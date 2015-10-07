#include "header.h"

int main(int argc, char *argv[])
{
  Parameters *params; // user defined parameters
  double ***phi;      // flux array
  FILE *fp = NULL;    // output file

  // Get inputs
  params = set_default_params();
  parse_params("parameters", params);
  print_params(params);

  // Open file for output
  if(params->print == 1){
    fp = fopen("output", "wb");
    fclose(fp);
  }

  // Initial guess of flux
  phi = init_flux(params);

  solve(phi, params);

  printf("keff = %f\n", params->k);

  // Free memory
  free_flux(phi);
  free(params);

  return 0;
}
