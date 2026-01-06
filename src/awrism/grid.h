#ifndef __RISM_GRID_H
#define __RISM_GRID_H

#include <stdlib.h>

struct Grid
{
  size_t n;
  double dr;
  double dk;
};
typedef struct Grid grid_t;

struct GridJac
{
  size_t n;
  double dr, dk, f, b, c;
};
typedef struct GridJac grid_jac_t;

struct GridEq
{
  size_t n;
  double dr, dk, f, b, c;
};
typedef struct GridEq grid_eq_t;

int ginit(int n, double dr, grid_t *g);
int grid_eq_mk(int n, double dr, grid_eq_t *p);
int grid_jac_mk(int n, float dr, grid_jac_t *p);

#endif /* __RISM_GRID_H */
