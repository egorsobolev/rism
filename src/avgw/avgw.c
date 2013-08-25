#include "avgw.h"

#include "shapes.h"
#include <rism/grid.h>
#include <rism/interp.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include <fftw3.h>

int avgw_func_init(const sgrid_t *g, avgw_func_t *f)
{
  f->np = g->np + 1;
  f->nz = 0;
  f->I = 0.0;
  f->s = (float *) calloc(3 * f->np, sizeof(float));
  if (!f->s)
    return -1;
  f->s2 = f->s + f->np;
  f->u = f->s2 + f->np;
  return 0;
}

void avgw_func_free(avgw_func_t *f)
{
  free(f->s);
}

int avgw_mtx_init(int n, gridparam_t *p, int nfun, avgw_mtx_t *W)
{
  int nb, i;
  nb = sizeof(float) + 2 * sizeof(int);
  for (i = 0; i < n; i++) {
    W[i].npcut = (int *) calloc(nfun, nb);
    if (!W[i].npcut)
      break;
    W[i].nz = W[i].npcut + nfun;
    W[i].Icut = (float *) (W[i].nz + nfun);
    W[i].nfun = nfun;
    W[i].dr = p[i].dr;
    W[i].np = p[i].np;
  }
  if (i == n)
    return 0;
  for (; i >= 0; i--)
    free(W[i].npcut);
}

void avgw_mtx_free(int n, avgw_mtx_t *W)
{
  int i;
  for (i = 0; i < n; i++)
    free(W[i].npcut);
}

int avgw_write_hdr(avgw_mtx_t *W, FILE *f)
{
  if (fwrite(&W->np, sizeof(int), 1, f) != 1)
    return -1;
  if (fwrite(&W->dr, sizeof(float), 1, f) != 1)
    return -1;
  if (fwrite(&W->nfun, sizeof(int), 1, f) != 1)
    return -1;
  if (fwrite(W->npcut, sizeof(int), W->nfun, f) != W->nfun)
    return -1;
  return 0;
}

int avgw_write_info(avgw_mtx_t *W, FILE *f)
{
  if (fwrite(W->nz, sizeof(int), W->nfun, f) != W->nfun)
    return -1;
  if (fwrite(W->Icut, sizeof(float), W->nfun, f) != W->nfun)
    return -1;
  return 0;
}

void avgw_hist2aw(sgrid_t *g, int n, int i0, int nsamp, const unsigned *h, avgw_func_t *f)
{
  float k, w0, w1;
  int m, i;
  /* expand histogram */
  memset(g->d, 0, (g->np - 1) * sizeof(float));
  m = i0 + n;
  if (m > g->np)
    m = g->np;
  k = 1.0 / nsamp;
  for (i = i0; i < m; i++)
    g->d[i - 1] = k * h[i - i0] / i;

  /* make fft */
  fftwf_execute(g->p);
 
  k = 0.5 * g->np / M_PI;
  f->s[0] = 1.0;
  for (i = 1; i < g->np; i++)
    f->s[i] = k * g->d[i - 1] / i;
  f->s[g->np] = 0.0;

  /* int |w(k)|dk, k=0..inf */
  f->I = 0.5; /* w(0) = 1 */ 
  f->nz = 0;
  w1 = 0.0;
  for (i = f->np - 2; i > 0; --i) {
    w0 = f->s[i];
    if (w1 == 0.0 || (w0 * w1) < 0.0) {
      f->I -= fabs(w0 * w1 / (w0 - w1));
      f->nz++;
    }
    f->I += fabs(w0);
    w1 = w0;
  }
}

void avgw_itail(const sgrid_t *g, const avgw_func_t *f, avgw_shapes_t *s, avgw_cutparam_t *c)
{
  int i, j, k, nz, nzmax;
  float w0, w1, lcut, lw;
  double I1, err;

  i = f->np - 2;
  nz = 0;
  w1 = 0.0;
  I1 = 0.0;

  nzmax = f->nz - AVGW_MINZEROS;

  for (k = 0; k < s->n; k++) {
    j = s->o[k];
    c[j].npcut = i;
    c[j].Icut = I1;
    c[j].nz = f->nz - nz;

    err = f->I * s->p[j].trun;
    /* int |w(k)|dk, k=icut..inf */
    while (i >= 0 && I1 < err && nz <= nzmax) {
      w0 = f->s[i];
      if (w1 == 0.0 || (w0 * w1) < 0.0) {
	nz++;
	c[j].Icut = I1;
	c[j].nz = f->nz - nz;
	I1 -= fabs(w0 * w1 / (w0 - w1));
	lcut = g->dt * (i + w0 / fabs(w0 - w1));
      }
      I1 += fabs(w0);
      w1 = w0;
      i--;
    }
    lw = (s->p[j].np - 1) * s->p[j].dt;
    /* test it */
    c[j].npcut = lcut < lw ? (int) (lcut / s->p[j].dt) - (fmodf(lcut, s->p[j].dt) <= FLT_EPSILON) : (s->p[j].np - 1);
  }
}

void avgw_reshape(const sgrid_t *g, const avgw_func_t *f, const gridparam_t *p, const avgw_cutparam_t *c, float *fcut)
{
  int i, k, each;
  float h;

  h = p->dt / g->dt;
  if (p->interp) {
    ssplint_uni(f->np, f->s, f->s2, c->npcut, h, h, fcut);
  } else {
    each = (int) (h + 0.5);
    k = each;
    for (i = 0; i < c->npcut; i++) {
      fcut[i] = f->s[k];
      k += each;
    }
  }
}
