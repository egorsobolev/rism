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
  g->d = (float *) fftwf_malloc((np - 1) * sizeof(float));
  if (!g->d)
    return -1;
  g->p = fftwf_plan_r2r_1d(np - 1, g->d, g->d, FFTW_RODFT00, FFTW_ESTIMATE);
  if (!g->p) {
    fftwf_free(g->d);
    return -1;
  }
  return 0;
}

void sgrid_free(sgrid_t *g)
{
  fftwf_destroy_plan(g->p);
  fftwf_free(g->d);
}

int dgrid_init(int np, double dr, dgrid_t *g)
{
  g->np = np;
  g->dr = dr;
  g->dt = M_PI / (dr * np);
  g->f = 4.0 * M_PI * dr;
  g->b = 2.0 / (g->f * np);
  g->d = (double *) fftw_malloc((np - 1) * sizeof(double));
  if (!g->d)
    return -1;
  g->p = fftw_plan_r2r_1d(np - 1, g->d, g->d, FFTW_RODFT00, FFTW_ESTIMATE);
  if (!g->p) {
    fftw_free(g->d);
    return -1;
  }
  return 0;
}

void dgrid_free(dgrid_t *g)
{
  fftw_destroy_plan(g->p);
  fftw_free(g->d);
}
