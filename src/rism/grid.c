#include <rism/grid.h>

#include <stdlib.h>

#define _USE_MATH_DEFINES
#include <math.h>

int grid_init(int n, double dr, grid_t *g)
{
  g->ngrid = n;
  g->n = n - 1;
  g->dr = dr;
  g->dk = M_PI / (dr * n);
  g->f = 2.0 * M_PI * dr * dr;
  g->b = 0.5 / (g->f * n);
  g->c = 0.5 / n;
  return 0;
}
