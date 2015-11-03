#include "header.h"

void parse_params(char *filename, Parameters *params)
{
  char line[256], *s;
  FILE *fp = fopen(filename, "r");

  while((s = fgets(line, sizeof(line), fp)) != NULL){

    if(line[0] == '#') continue;
    s = strtok(line, "=");
    if(s == NULL) continue;

    // Number of grid points (n_grid)
    else if(strcmp(s, "n_grid") == 0)
      params->n_grid = atoi(strtok(NULL, "=\n"));

    // Grid spacing (h)
    else if(strcmp(s, "h") == 0)
      params->h = atof(strtok(NULL, "=\n"));

    // Number of energy groups (groups)
    else if(strcmp(s, "groups") == 0)
      params->G = atoi(strtok(NULL, "=\n"));

    // Average cos of scattering angle (mu)
    else if(strcmp(s, "mu") == 0)
      params->mu = atof(strtok(NULL, "=\n"));

    // Boundary conditions (bc)
    else if(strcmp(s, "bc") == 0){
      s = strtok(NULL, "=\n");
      if(strcasecmp(s, "vacuum") == 0)
        params->bc = VACUUM;
      else if(strcasecmp(s, "periodic") == 0)
        params->bc = PERIODIC;
      else print_error("Invalid boundary condition");
    }

    // Maximum number of inner iterations (max_inner)
    else if(strcmp(s, "max_inner") == 0)
      params->max_inner = atoi(strtok(NULL, "=\n"));

    // Maximum number of outer iterations (max_outer)
    else if(strcmp(s, "max_outer") == 0)
      params->max_outer = atoi(strtok(NULL, "=\n"));

    // Stopping condition (thresh)
    else if(strcmp(s, "thresh") == 0)
      params->thresh = atof(strtok(NULL, "=\n"));

    // Whether to output flux (write_flux)
    else if(strcmp(s, "write_flux") == 0){
      s = strtok(NULL, "=\n");
      if(strcasecmp(s, "true") == 0)
        params->write_flux = TRUE;
      else if(strcasecmp(s, "false") == 0)
        params->write_flux = FALSE;
      else
        print_error("Invalid option for parameter 'write_flux': must be 'true' or 'false'");
    }

    // Path to write flux to (flux_file)
    else if(strcmp(s, "flux_file") == 0){
      s = strtok(NULL, "=\n");
      params->flux_file = malloc(strlen(s)*sizeof(char)+1);
      strcpy(params->flux_file, s);
    }

    // Path to read group constants from (group_file)
    else if(strcmp(s, "group_file") == 0){
      s = strtok(NULL, "=\n");
      params->group_file = malloc(strlen(s)*sizeof(char)+1);
      strcpy(params->group_file, s);
    }

    // Unknown option
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

    // Number of energy groups (-e)
    else if(strcmp(arg, "-e") == 0){
      if(++i < argc) params->G = atoi(argv[i]);
      else print_error("Error reading command line input '-e'");
    }

    // Average cos of scattering angle (-m)
    else if(strcmp(arg, "-m") == 0){
      if(++i < argc) params->mu = atof(argv[i]);
      else print_error("Error reading command line input '-m'");
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

    // Path to read group constants from (-q)
    else if(strcmp(arg, "-q") == 0){
      if(++i < argc){
        if(params->group_file != NULL) free(params->group_file);
        params->group_file = malloc(strlen(argv[i])*sizeof(char)+1);
        strcpy(params->group_file, argv[i]);
      }
      else print_error("Error reading command line input '-p'");
    }

    // Unknown option
    else print_error("Error reading command line input");
  }

  // Read in group constants from separate file
  if(params->group_file == NULL)
    print_error("Must specify file to read in group constants");
  set_group_constants(params);

  // Set remaining parameters
  params->L = params->h*params->n_grid;
  params->k = 1;

  // Validate inputs
  if(params->write_flux == TRUE && params->flux_file == NULL)
    params->flux_file = "flux.dat";
  if(params->n_grid <= 0)
    print_error("Number of grid points must be greater than zero");
  if(params->h <= 0)
    print_error("Grid spacing must be greater than zero");
  if(params->G <= 0)
    print_error("Number of energy groups must be greater than zero");
  if(params->mu < -1 || params->mu > 1)
    print_error("mu must be in range [-1, 1]");
  if(params->max_inner <= 0 || params->max_outer <= 0)
    print_error("Maximum number of iterations must be greater than zero");
  if(params->thresh <= 0)
    print_error("Threshold must be greater than zero");

  return;
}

void set_group_constants(Parameters *params)
{
  char line[512], *s;
  FILE *fp = fopen(params->group_file, "r");
  int i, j;
  double *arr;
  double *xs_s;

  // Allocate memory for group constants
  params->xs_f = calloc(params->G, sizeof(double));
  params->xs_a = calloc(params->G, sizeof(double));
  params->xs_s = malloc(params->G*sizeof(double*));
  arr = calloc(params->G*params->G, sizeof(double));
  for(i=0; i<params->G; i++){
    params->xs_s[i] = arr;
    arr += params->G;
  }
  params->xs_t = calloc(params->G, sizeof(double));
  params->xs_r = calloc(params->G, sizeof(double));
  params->nu = calloc(params->G, sizeof(double));
  params->D = calloc(params->G, sizeof(double));
  params->chi = calloc(params->G, sizeof(double));

  while((s = fgets(line, sizeof(line), fp)) != NULL){

    s = strtok(line, "\n");

    // Average number of fission neutrons produced
    if(strcmp(s, "nu") == 0){
      s = fgets(line, sizeof(line), fp);
      s = strtok(line, " \n");
      for(i=0; i<params->G; i++){
        params->nu[i] = atof(s);
        s = strtok(NULL, " \n");
      }
    }

    // Probability of fission neutron being born in group
    else if(strcmp(s, "chi") == 0){
      s = fgets(line, sizeof(line), fp);
      s = strtok(line, " \n");
      for(i=0; i<params->G; i++){
        params->chi[i] = atof(s);
        s = strtok(NULL, " \n");
      }
    }

    // Fission cross section
    else if(strcmp(s, "xs_f") == 0){
      s = fgets(line, sizeof(line), fp);
      s = strtok(line, " \n");
      for(i=0; i<params->G; i++){
        params->xs_f[i] = atof(s);
        s = strtok(NULL, " \n");
      }
    }

    // Absorption cross section
    else if(strcmp(s, "xs_a") == 0){
      s = fgets(line, sizeof(line), fp);
      s = strtok(line, " \n");
      for(i=0; i<params->G; i++){
        params->xs_a[i] = atof(s);
        s = strtok(NULL, " \n");
      }
    }
    
    // Scattering cross section
    else if(strcmp(s, "xs_s") == 0){
      for(i=0; i<params->G; i++){
        s = fgets(line, sizeof(line), fp);
        s = strtok(line, " \n");
        for(j=0; j<params->G; j++){
          params->xs_s[i][j] = atof(s);
          s = strtok(NULL, " \n");
        }
      }
    }

    else printf("Unknown value '%s' in group constants file.\n", s);
  }

  // Set remaining group constants
  xs_s = calloc(params->G, sizeof(double));
  for(i=0; i<params->G; i++){
    for(j=0; j<params->G; j++){
      xs_s[i] += params->xs_s[i][j];
    }
    params->xs_t[i] = params->xs_f[i] + params->xs_a[i] + xs_s[i];
    params->xs_r[i] = params->xs_t[i] - params->xs_s[i][i];
    params->D[i] = 1/(3*params->xs_t[i] - params->mu*xs_s[i]);
  }
  free(xs_s);

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
  printf("Number of energy groups:              %d\n", params->G);
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
void write_flux(double ****phi, Parameters *params, FILE *fp)
{
  int i, j, k;
  double ***S;
  double S_max;

  S = matrix3D(params->n_grid, params->n_grid, params->n_grid);
  compute_source(phi, S, params, &S_max);
  fp = fopen(params->flux_file, "w");

  for(i=0; i<params->n_grid; i++){
    for(j=0; j<params->n_grid; j++){
      for(k=0; k<params->n_grid; k++){
        fprintf(fp, "%e ", S[i][j][k]);
      }
      fprintf(fp, "\n");
    }
  }

  fclose(fp);
  free_matrix3D(S);

  return;
}
