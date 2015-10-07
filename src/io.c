#include "header.h"

void print_slice(double ***phi, Parameters *params, FILE *fp)
{
  int i, j, k;

  fp = fopen("output", "a");
  i = params->n_grid/2;
  for(j=1; j<params->n_grid+1; j++){
    for(k=1; k<params->n_grid+1; k++){
      fprintf(fp, "%f ", phi[i][j][k]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
}

void parse_params(char *filename, Parameters *params)
{
  char line[256], *s;
  FILE *fp = fopen(filename, "r");

  while((s = fgets(line, sizeof(line), fp)) != NULL){

    if(line[0] == '#') continue;
    s = strtok(line, "=");
    if(s == NULL) continue;

    // Set parameters
    else if(strcmp(s, "n_grid") == 0) params->n_grid = atoi(strtok(NULL, "=\n"));
    else if(strcmp(s, "h") == 0) params->h = atof(strtok(NULL, "=\n"));
    else if(strcmp(s, "macro_xs_f") == 0) params->macro_xs_f = atof(strtok(NULL, "=\n"));
    else if(strcmp(s, "macro_xs_a") == 0) params->macro_xs_a = atof(strtok(NULL, "=\n"));
    else if(strcmp(s, "macro_xs_e") == 0) params->macro_xs_e = atof(strtok(NULL, "=\n"));
    else if(strcmp(s, "mu") == 0) params->mu = atof(strtok(NULL, "=\n"));
    else if(strcmp(s, "nu") == 0) params->nu = atof(strtok(NULL, "=\n"));
    else if(strcmp(s, "bc") == 0){
      s = strtok(NULL, "=\n");
      if(strcasecmp(s, "vacuum") == 0) params->bc = 0;
      else if(strcasecmp(s, "reflective") == 0) params->bc = 1;
      else if(strcasecmp(s, "periodic") == 0) params->bc = 2;
      else print_error("Invalid boundary condition");
    }
    else if(strcmp(s, "max_inner") == 0) params->max_inner = atoi(strtok(NULL, "=\n"));
    else if(strcmp(s, "max_outer") == 0) params->max_outer = atoi(strtok(NULL, "=\n"));
    else if(strcmp(s, "thresh") == 0) params->thresh = atof(strtok(NULL, "=\n"));
    else if(strcmp(s, "print") == 0) params->print = atoi(strtok(NULL, "=\n"));
    else printf("Unknown value '%s' in config file.\n", s);
  }

  // Set remaining parameters
  params->L = params->h*params->n_grid;
  params->macro_xs_t = params->macro_xs_f + params->macro_xs_a + params->macro_xs_e;
  params->D = 1/(3*params->macro_xs_t - params->mu*params->macro_xs_e);
  params->k = 1;

  // Validate inputs
  if(params->n_grid <= 0)
    print_error("Number of grid points must be greater than zero");
  if(params->h <= 0)
    print_error("Grid spacing must be greater than zero");
  if(params->macro_xs_f < 0 || params->macro_xs_a < 0 || params->macro_xs_e < 0)
    print_error("Macroscopic cross section values cannot be negative");
  if(params->mu < -1 || params->mu > 1)
    print_error("mu must be in range [-1, 1]");
  if(params->nu < 0)
    print_error("nu cannot be negative");
  if(params->max_inner <= 0 || params->max_outer <= 0)
    print_error("Maximum number of iterations must be greater than zero");
  if(params->thresh <= 0)
    print_error("Threshold must be greater than zero");
}

void print_params(Parameters *params)
{
  char *bc = NULL;
  if(params->bc == 0) bc = "Vacuum";
  else if(params->bc == 1) bc = "Reflective";
  else if(params->bc == 2) bc = "Periodic";
  border_print();
  center_print("INPUT SUMMARY", 79);
  border_print();
  printf("Number of grid points:                %d\n", params->n_grid);
  printf("Grid spacing:                         %f\n", params->h);
  printf("Detector length:                      %f\n", params->L);
  printf("Average cos of scattering angle:      %f\n", params->mu);
  printf("Average number of fission neutrons:   %f\n", params->nu);
  printf("Diffusion coefficient:                %f\n", params->D);
  printf("Maximum number of inner iterations:   %d\n", params->max_inner);
  printf("Maximum number of outer iterations:   %d\n", params->max_outer);
  printf("Stopping threshold:                   %f\n", params->thresh);
  printf("Boundary conditions:                  %s\n", bc);
  border_print();
}

void print_error(char *message) 
{ 
  printf("ERROR: %s\n", message); 
  exit(1); 
}

void border_print(void)
{
  printf( "=========================================="
     "======================================\n");
}

// Prints comma separated integers - for ease of reading
void fancy_int(long a)
{
  if(a < 1000)
    printf("%ld\n",a);
  else if(a >= 1000 && a < 1000000)
    printf("%ld,%03ld\n", a/1000, a % 1000);
  else if(a >= 1000000 && a < 1000000000)
    printf("%ld,%03ld,%03ld\n", a/1000000, (a % 1000000)/1000, a % 1000);
  else if(a >= 1000000000)
    printf("%ld,%03ld,%03ld,%03ld\n", a / 1000000000,
       (a % 1000000000)/1000000, (a % 1000000)/1000, a % 1000);
  else
    printf("%ld\n",a);
}

// Prints Section titles in center of 80 char terminal
void center_print(const char *s, int width)
{
  int length = strlen(s);
  int i;
  for (i=0; i<=(width-length)/2; i++){
    fputs(" ", stdout);
  }
  fputs(s, stdout);
  fputs("\n", stdout);
}
