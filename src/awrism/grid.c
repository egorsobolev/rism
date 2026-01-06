#include "grid.h"

#include <stdlib.h>

#define _USE_MATH_DEFINES
#include <math.h>

int ginit(int n, double dr, grid_t *g)
{
  g->n = n;
  g->dr = dr;
  g->dk = M_PI / (dr * n);
  return 0;
}

int grid_eq_mk(int n, double dr, grid_eq_t *p)
{
  p->n = n;
  n += 1;
  p->dr = dr;
  p->dk = M_PI / (dr * n);
  p->f = 2.0 * M_PI * dr * dr;
  p->b = 0.5 / (p->f * n);
  p->c = 0.5 / n;
  return 0;
}

int grid_jac_mk(int n, float dr, grid_jac_t *p)
{
  p->n = n;
  n += 1;
  p->dr = dr;
  p->dk = M_PI / (dr * n);
  p->f = 2.0 * M_PI * dr * dr;
  p->b = 0.5 / (p->f * n);
  p->c = 0.5 / n;
  return 0;
}
