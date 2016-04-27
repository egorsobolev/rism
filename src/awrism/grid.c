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
  p->data = (double *) fftw_malloc(sizeof(double) * n);
  if (!p->data)
    return -1;
  p->plan = fftw_plan_r2r_1d(n, p->data, p->data, FFTW_RODFT00, FFTW_ESTIMATE);
  if (!p->plan) {
    fftw_free(p->data);
    return -2;
  }
  p->n = n;
  n += 1;
  p->dr = dr;
  p->dk = M_PI / (dr * n);
  p->f = 2.0 * M_PI * dr * dr;
  p->b = 0.5 / (p->f * n);
  p->c = 0.5 / n;
  return 0;
}

void grid_eq_del(grid_eq_t *p)
{
  fftw_destroy_plan(p->plan);
  fftw_free(p->data);
}

int grid_jac_mk(int n, float dr, grid_jac_t *p)
{
  p->data = (float *) fftwf_malloc(sizeof(float) * n);
  if (!p->data)
    return -1;
  p->plan = fftwf_plan_r2r_1d(n, p->data, p->data, FFTW_RODFT00, FFTW_ESTIMATE);
  if (!p->plan) {
    fftwf_free(p->data);
    return -2;
  }
  p->n = n;
  n += 1;
  p->dr = dr;
  p->dk = M_PI / (dr * n);
  p->f = 2.0 * M_PI * dr * dr;
  p->b = 0.5 / (p->f * n);
  p->c = 0.5 / n;
  return 0;
}

void grid_jac_del(grid_jac_t *p)
{
  fftwf_destroy_plan(p->plan);
  fftwf_free(p->data);
}
