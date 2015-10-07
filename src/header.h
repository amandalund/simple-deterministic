#ifndef HEADER
#define HEADER

#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<time.h>
#include<sys/time.h>
#include<math.h>
#include<string.h>

#define FLAG -13

// Boundary conditions
#define VACUUM 0
#define PERIODIC 1

// User defined parameters
typedef struct Parameters_{
  int n_grid;                // number of grid points in x, y, and z
  double L;                  // detector length
  double h;                  // grid spacing
  double macro_xs_a;         // absorption macro xs
  double macro_xs_e;         // elastic scattering  macro xs
  double macro_xs_f;         // fission macro xs
  double macro_xs_t;         // total macro xs
  double mu;                 // average cos of scattering angle in scattering collision
  double D;                  // diffusion coefficient
  double nu;                 // average number of fission neutrons produced
  double k;                  // multiplication factor
  int bc;                    // boundary conditions
  int max_inner;             // maximum number of inner iterations
  int max_outer;             // maximum number of outer iterations
  double thresh;             // threshold for termination condition
  int print;                 // whether to print output
} Parameters;

// initialize.c function prototypes
Parameters *set_default_params(void);
double ***init_flux(Parameters *params);
void free_flux(double ***phi);

// utils.c function prototypes
double timer(void);
double ***matrix3D(size_t l, size_t m, size_t n);
void free_matrix3D(double ***m);

// io.c function prototypes
void print_slice(double ***phi, Parameters *params, FILE *fp);
void parse_params(char *filename, Parameters *params);
void print_params(Parameters *params);
void print_error(char *message);
void border_print(void);
void fancy_int(long a);
void center_print(const char *s, int width);

// solvers.c function prototypes
void solve(double ***phi0, Parameters *params);

#endif