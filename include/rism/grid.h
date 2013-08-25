#ifndef __RISM_GRID_H
#define __RISM_GRID_H

#include <fftw3.h>

struct GRID_FLOAT
{
  int np;
  float dr;
  float dt;
  float f;
  float b;
  fftwf_plan p;
  float *d;
};

typedef struct GRID_FLOAT sgrid_t;

struct GRID_DOUBLE
{
  int np;
  double dr;
  double dt;
  double f;
  double b;
  fftw_plan p;
  double *d;
};

typedef struct GRID_DOUBLE dgrid_t;

int sgrid_init(int np, float dr, sgrid_t *g);
int dgrid_init(int np, double dr, dgrid_t *g);
void sgrid_free(sgrid_t *g);
void dgrid_free(dgrid_t *g);

#endif /* __RISM_GRID_H */
