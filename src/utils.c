#include "header.h"

double timer(void)
{
  struct timeval time;
  gettimeofday(&time, 0);
  return time.tv_sec + time.tv_usec/1000000.0;
}

double ***matrix3D(size_t l, size_t m, size_t n)
{
  int i, j;
  double *c = malloc(l*m*n*sizeof(double));
  double **b = malloc(l*m*sizeof(double*));
  double ***a = malloc(l*sizeof(double**));

  for(i=0; i<l; i++){
    a[i] = b;
    b += m;
    for(j=0; j<m; j++){
      a[i][j] = c;
      c += n;
    }
  }

  return a;
}

void free_matrix3D(double ***m)
{
  free(m[0][0]);
  free(m[0]);
  free(m);

  return;
}
