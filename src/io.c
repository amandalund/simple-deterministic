#include "header.h"

void parse_params(char *filename, Parameters *params)
{
  char line[256], *s;
  FILE *fp = fopen(filename, "r");

  while((s = fgets(line, sizeof(line), fp)) != NULL){

    if(line[0] == '#') continue;
    s = strtok(line, "=");
    if(s == NULL) continue;

    // Set parameters
    else if(strcmp(s, "n_grid") == 0)
      params->n_grid = atoi(strtok(NULL, "=\n"));
    else if(strcmp(s, "h") == 0)
      params->h = atof(strtok(NULL, "=\n"));
    else if(strcmp(s, "macro_xs_f") == 0)
      params->macro_xs_f = atof(strtok(NULL, "=\n"));
    else if(strcmp(s, "macro_xs_a") == 0)
      params->macro_xs_a = atof(strtok(NULL, "=\n"));
    else if(strcmp(s, "macro_xs_e") == 0)
      params->macro_xs_e = atof(strtok(NULL, "=\n"));
    else if(strcmp(s, "mu") == 0)
      params->mu = atof(strtok(NULL, "=\n"));
    else if(strcmp(s, "nu") == 0)
      params->nu = atof(strtok(NULL, "=\n"));
    else if(strcmp(s, "bc") == 0){
      s = strtok(NULL, "=\n");
      if(strcasecmp(s, "vacuum") == 0)
        params->bc = VACUUM;
      else if(strcasecmp(s, "periodic") == 0)
        params->bc = PERIODIC;
      else print_error("Invalid boundary condition");
    }
    else if(strcmp(s, "max_inner") == 0)
      params->max_inner = atoi(strtok(NULL, "=\n"));
    else if(strcmp(s, "max_outer") == 0)
      params->max_outer = atoi(strtok(NULL, "=\n"));
    else if(strcmp(s, "thresh") == 0)
      params->thresh = atof(strtok(NULL, "=\n"));
    else if(strcmp(s, "write_flux") == 0){
      s = strtok(NULL, "=\n");
      if(strcasecmp(s, "true") == 0)
        params->write_flux = TRUE;
      else if(strcasecmp(s, "false") == 0)
        params->write_flux = FALSE;
      else
        print_error("Invalid option for parameter 'write_flux': must be 'true' or 'false'");
    }
    else if(strcmp(s, "flux_file") == 0){
      s = strtok(NULL, "=\n");
      params->flux_file = malloc(strlen(s)*sizeof(char)+1);
      strcpy(params->flux_file, s);
    }
    else printf("Unknown value '%s' in config file.\n", s);
  }

  return;
}

void read_CLI(int argc, char *argv[], Parameters *params)
{
  int i;
  char *arg;

  // Collect raw input
  for(i=1; i<argc; i++){
    arg = argv[i];

    // Number of grid points (-g)
    if(strcmp(arg, "-g") == 0){
      if(++i < argc) params->n_grid = atoi(argv[i]);
      else print_error("Error reading command line input '-g'");
    }

    // Boundary conditions (-c)
    else if(strcmp(arg, "-c") == 0){
      if(++i < argc){
        if(strcasecmp(argv[i], "vacuum") == 0)
          params->bc = VACUUM;
        else if(strcasecmp(argv[i], "periodic") == 0)
          params->bc = PERIODIC;
        else
          print_error("Invalid boundary condition");
      }
      else print_error("Error reading command line input '-c'");
    }

    // Grid spacing (-h)
    else if(strcmp(arg, "-h") == 0){
      if(++i < argc) params->h = atof(argv[i]);
      else print_error("Error reading command line input '-h'");
    }

    // Absorption macro xs (-a)
    else if(strcmp(arg, "-a") == 0){
      if(++i < argc) params->macro_xs_a = atof(argv[i]);
      else print_error("Error reading command line input '-a'");
    }

    // Elastic macro xs (-e)
    else if(strcmp(arg, "-e") == 0){
      if(++i < argc) params->macro_xs_e = atof(argv[i]);
      else print_error("Error reading command line input '-e'");
    }

    // Fission macro xs (-f)
    else if(strcmp(arg, "-f") == 0){
      if(++i < argc) params->macro_xs_f = atof(argv[i]);
      else print_error("Error reading command line input '-f'");
    }

    // Average cos of scattering angle (-m)
    else if(strcmp(arg, "-m") == 0){
      if(++i < argc) params->mu = atof(argv[i]);
      else print_error("Error reading command line input '-m'");
    }

    // Average number of fission neutrons produced (-n)
    else if(strcmp(arg, "-n") == 0){
      if(++i < argc) params->nu = atof(argv[i]);
      else print_error("Error reading command line input '-n'");
    }

    // Maximum number of inner iterations (-i)
    else if(strcmp(arg, "-i") == 0){
      if(++i < argc) params->max_inner = atoi(argv[i]);
      else print_error("Error reading command line input '-i'");
    }

    // Maximum number of outer iterations (-o)
    else if(strcmp(arg, "-o") == 0){
      if(++i < argc) params->max_outer = atoi(argv[i]);
      else print_error("Error reading command line input '-o'");
    }

    // Stopping condition (-t)
    else if(strcmp(arg, "-t") == 0){
      if(++i < argc) params->thresh = atof(argv[i]);
      else print_error("Error reading command line input '-t'");
    }

    // Whether to output flux (-w)
    else if(strcmp(arg, "-w") == 0){
      if(++i < argc){
        if(strcasecmp(argv[i], "true") == 0)
          params->write_flux = TRUE;
        else if(strcasecmp(argv[i], "false") == 0)
          params->write_flux = FALSE;
        else
          print_error("Invalid option for parameter 'write_flux': must be 'true' or 'false'");
      }
      else print_error("Error reading command line input '-w'");
    }

    // Path to write flux to (-p)
    else if(strcmp(arg, "-p") == 0){
      if(++i < argc){
        if(params->flux_file != NULL) free(params->flux_file);
        params->flux_file = malloc(strlen(argv[i])*sizeof(char)+1);
        strcpy(params->flux_file, argv[i]);
      }
      else print_error("Error reading command line input '-p'");
    }

    else print_error("Error reading command line input");
  }

  // Set remaining parameters
  params->L = params->h*params->n_grid;
  params->macro_xs_t = params->macro_xs_f + params->macro_xs_a + params->macro_xs_e;
  params->D = 1/(3*params->macro_xs_t - params->mu*params->macro_xs_e);
  params->k = 1;

  // Validate inputs
  if(params->write_flux == TRUE && params->flux_file == NULL)
    params->flux_file = "flux.dat";
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

  return;
}

void print_params(Parameters *params)
{
  char *bc = NULL;
  if(params->bc == VACUUM) bc = "Vacuum";
  else if(params->bc == PERIODIC) bc = "Periodic";
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
  printf("Stopping threshold:                   %e\n", params->thresh);
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

// Prints solution to file
void write_flux(double ***phi, Parameters *params, FILE *fp)
{
  int i, j, k;

  fp = fopen(params->flux_file, "w");

  for(i=0; i<params->n_grid; i++){
    for(j=0; j<params->n_grid; j++){
      for(k=0; k<params->n_grid; k++){
        fprintf(fp, "%e ", phi[i][j][k]);
      }
      fprintf(fp, "\n");
    }
  }

  fclose(fp);

  return;
}
