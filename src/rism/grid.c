#include <rism/grid.h>

#define _USE_MATH_DEFINES
#include <math.h>

int sgrid_init(int np, float dr, sgrid_t *g)
{
  g->np = np;
  g->dr = dr;
  g->dt = M_PI / (dr * np);
  g->f = 4.0 * M_PI * dr;
  g->b = 2.0 / (g->f * np);
  return 0;
}

int dgrid_init(int np, double dr, dgrid_t *g)
{
  g->np = np;
  g->dr = dr;
  g->dt = M_PI / (dr * np);
  g->f = 4.0 * M_PI * dr;
  g->b = 2.0 / (g->f * np);
  return 0;
}
