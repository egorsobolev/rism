#ifndef __RISM_GRID_H
#define __RISM_GRID_H

#include <fftw3.h>

struct Grid
{
  int n;
  double dr;
  double dk;
};
typedef struct Grid grid_t;

struct GridJac
{
  fftwf_plan plan;
  float *data;
  int n;
  double dr, dk, f, b, c;
};
typedef struct GridJac grid_jac_t;

struct GridEq
{
  fftw_plan plan;
  double *data;
  int n;
  double dr, dk, f, b, c;
};
typedef struct GridEq grid_eq_t;

int ginit(int n, double dr, grid_t *g);
int grid_eq_mk(int n, double dr, grid_eq_t *p);
void grid_eq_del(grid_eq_t *p);
int grid_jac_mk(int n, float dr, grid_jac_t *p);
void grid_jac_del(grid_jac_t *p);

#endif /* __RISM_GRID_H */
