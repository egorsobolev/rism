#include "eqoz.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <cblas.h>

void rismaw_getx(const void *prm, const double *d, float *x)
{
  int i, j, l, k;
  const rismaw_t *p = (rismaw_t *) prm;

  l = 0;
  for (i = 0; i < p->nfun; i++) {
    k = i * p->ge.n;
    for (j = 0; j < p->gj.n; j++) {
      x[l] = (float) d[k + j];
      l++;
    }
  }
}

void rismaw_putx(const void *prm, const float *x, double *d)
{
  int i, j, l, k;
  const rismaw_t *p = (rismaw_t *) prm;

  l = 0;
  for (i = 0; i < p->nfun; i++) {
    k = i * p->ge.n;
    for (j = 0; j < p->gj.n; j++) {
      d[k + j] = (double) x[l];
      l++;
    }
  }
}

#define max(a, b) ((a) > (b) ? a : b)
#define min(a, b) ((a) < (b) ? a : b)

int rismaw_eq(void *prm, const double *tuv, double *d, double *en)
{
  double *r;
  int i, j, u, v, k, l, s, np, incu, i0, n;
  double a;
  
  rismaw_t *p = (rismaw_t *) prm;
  grid_eq_t *g = &p->ge;

  incu = p->natv * g->n;
  np = p->natu * incu;

  r = (double *) malloc(np * sizeof(double));
  if (!r) {
    return -1;
  }

  /* cuv = C[tuv(r)] */
  p->closure(p, tuv, r, p->dcdt, en);

  /* r = fft(cuv) */
  l = 0;
  n = min(p->puv.ngrid, g->n);
  for (u = 0; u < p->natu; u++) {
    k = p->puv.atyp[u] * p->natv * p->puv.ngrid;
    for (v = 0; v < p->natv; v++) {
      for (i = 0; i < n; i++)
        g->data[i] = (r[l + i] - p->puv.asympr[k + i]) * (i + 1);
      for (; i < g->n; i++)
        g->data[i] = 0.0;

      fftw_execute(g->plan);

      for (i = 0; i < n; i++)
        r[l + i] = (g->f * g->data[i] + p->puv.asympk[k + i]) * p->v.symc[v];
      for (; i < g->n; i++)
        r[l + i] = 0.0;

      l += g->n;
      k += p->puv.ngrid;
    }
  }

  /* d = wuu * r = wuu * fft(cuv) */
  avgw_mul_dbl(&p->wuu, g->n, p->natv, p->natu, r, d);

  /* r = d * xvv - r = wuu * fft(cuv) * xvv - fft(cuv) */
  cblas_dscal(np, -1.0, r, 1);
  n = min(p->v.ngrid - 1, g->n);
  i0 = 1;
  for (v = 0; v < p->natv; v++) {
    /* side elements */
    for (j = 0; j < v; j++) {
      l = v * g->n;
      k = j * g->n;
      for (u = 0; u < p->natu; u++) {
        for (i = 0; i < n; i++) {
          a = p->v.xvv[i0 + i];
          r[l + i] += a * d[k + i];
          r[k + i] += a * d[l + i];
        }
        l += incu;
        k += incu;
      }
      i0 += p->v.ngrid + 1;
    }
    /* diagonal */
    l = v * g->n;
    for (u = 0; u < p->natu; u++) {
      for (i = 0; i < n; i++)
        r[l + i] += d[l + i] * p->v.xvv[i0 + i];
      l += incu;
    }
    i0 += p->v.ngrid + 1;
  }

  /* d = ifft(r) = ifft(wuu * fft(cuv) * xvv - fft(cuv)) */
  n = min(p->puv.ngrid, g->n);
  l = 0;
  for (u = 0; u < p->natu; u++) {
    k = p->puv.atyp[u] * p->natv * p->puv.ngrid;
    for (v = 0; v < p->natv; v++) {
      for (i = 0; i < n; i++)
        g->data[i] = r[l + i] / p->v.symc[v] + p->puv.asympk[k + i];
      for (; i < g->n; i++)
        g->data[i] = 0.0;

      fftw_execute(g->plan);

      for (i = 0; i < n; i++)
        d[l + i] = g->b * g->data[i] / (i + 1) - p->puv.asympr[k + i];
      for (; i < g->n; i++)
        d[l + i] = 0.0;

      l += g->n;
      k += p->puv.ngrid;
    }
  }

  /* d = d - tuv = ifft(wuu * fft(cuv) * xvv - fft(cuv)) - tuv */
  cblas_daxpy(np, -1.0, tuv, 1, d, 1);

  free(r);

  return 0;
}

int rismaw_Jx(void *prm, const float *x, float *r)
{
  float *d;
  int i, j, u, v, k, l, s, np, incu, i0, n;
  float a;

  rismaw_t *p = (rismaw_t *) prm;
  grid_jac_t *g = &p->gj;

  incu = p->natv * g->n;
  np = p->natu * incu;

  d = (float *) malloc(np * sizeof(float));
  if (!d) {
    return -1;
  }

  /* Jx = fft(w * fft(dcdt * x) * xvv - fft(dcdt * x)) - x */

  /* r1 = fft(dcdt * x) */
  l = 0;
  for (u = 0; u < p->natu; u++) {
    for (v = 0; v < p->natv; v++) {
      for (i = 0; i < g->n; i++)
        g->data[i] = p->dcdt[l + i] * x[l + i] * (i + 1);

      fftwf_execute(g->plan);

      a = g->f * (float) p->v.symc[v];
      for (i = 0; i < g->n; i++)
        r[l + i] = a * g->data[i];

      l += g->n;
    }
  }

  /* d = wuu * r = wuu * fft(cuv) */
  avgw_mul_flt(&p->wuu, p->reduc, g->n, p->natv, p->natu, r, d);

  /* r = d * xvv - r = wuu * fft(dcdt * x) * xvv - fft(dcdt * x) */
  cblas_sscal(np, -1.0f, r, 1);
  n = min(p->v.ngrid / p->reduc - 1, g->n);
  i0 = 1;
  for (v = 0; v < p->natv; v++) {
    /* side elements */
    for (j = 0; j < v; j++) {
      l = v * g->n;
      k = j * g->n;
      for (u = 0; u < p->natu; u++) {
        for (i = 0; i < n; i++) {
          a = (float) p->v.xvv[i0 + (i + 1) * p->reduc - 1];
          r[l + i] += a * d[k + i];
          r[k + i] += a * d[l + i];
        }
        l += incu;
        k += incu;
      }
      i0 += p->v.ngrid + 1;
    }
    /* diagonal */
    l = v * g->n;
    for (u = 0; u < p->natu; u++) {
      for (i = 0; i < n; i++)
        r[l + i] += d[l + i] * (float) p->v.xvv[i0 + (i + 1) * p->reduc - 1];
      l += incu;
    }
    i0 += p->v.ngrid + 1;
  }


  /* r = ifft(r3) - x = ifft(wuu * fft(dcdt * x) * xvv - fft(dcdt * x)) - x */
  l = 0;
  for (u = 0; u < p->natu; u++) {
    for (v = 0; v < p->natv; v++) {
      for (i = 0; i < g->n; i++)
        g->data[i] = r[l + i] / (float) p->v.symc[v];

      fftwf_execute(g->plan);

      for (i = 0; i < g->n; i++)
        r[l + i] = g->b * g->data[i] / (i + 1) - x[l + i];

      l += g->n;
    }
  }

  free(d);

  return 0;
}
