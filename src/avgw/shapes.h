#ifndef __RISM_AVGW_SHAPES_H
#define __RISM_AVGW_SHAPES_H

#include <rism/grid.h>
#include <stdio.h>

#ifndef GRID_SHAPES_MAX
#define GRID_SHAPES_MAX 16
#endif

struct GRIDPARAM
{
  int np;
  double dr;
  double dt;
  double trun;
  double l;
  int interp;
};
typedef struct GRIDPARAM gridparam_t;

struct AVGW_SHAPES
{
  int n;
  int interp;
  gridparam_t p[GRID_SHAPES_MAX];
  int o[GRID_SHAPES_MAX];
  FILE *f[GRID_SHAPES_MAX];
};
typedef struct AVGW_SHAPES avgw_shapes_t;

int avgw_shape_parse(const sgrid_t *g, const char *shape, gridparam_t *p);
void avgw_shapes_print(avgw_shapes_t *s);
int avgw_shapes_open(const char *prefix, int n, const gridparam_t *p, FILE **f);
int avgw_shapes_close(int n, FILE **f);

#endif /* __RISM_AVGW_SHAPES_H */
