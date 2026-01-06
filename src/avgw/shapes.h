#ifndef __RISM_AVGW_SHAPES_H
#define __RISM_AVGW_SHAPES_H

#include <rism/grid.h>
#include <stdio.h>
#ifdef MPI
#include <mpi.h>
#endif

#ifndef GRID_SHAPES_MAX
#define GRID_SHAPES_MAX 16
#endif

struct GRIDPARAM
{
  int np;
  double dr;
  double dk;
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
#ifdef MPI
  MPI_File f[GRID_SHAPES_MAX];
#else
  FILE *f[GRID_SHAPES_MAX];
#endif
};
typedef struct AVGW_SHAPES avgw_shapes_t;

int avgw_shape_parse(const grid_t *g, const char *shape, gridparam_t *p);
void avgw_shapes_print(avgw_shapes_t *s);
#ifdef MPI
int avgw_shapes_open(const char *prefix, int n, const gridparam_t *p, MPI_File *f);
int avgw_shapes_close(int n, MPI_File *f);
#else
int avgw_shapes_open(const char *prefix, int n, const gridparam_t *p, FILE **f);
int avgw_shapes_close(int n, FILE **f);
#endif

#endif /* __RISM_AVGW_SHAPES_H */
