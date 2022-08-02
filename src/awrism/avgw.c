#include "avgw.h" 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>

int avgw_read(const char *fn, avgw_t *w, int nge, int reduc)
{
  FILE *f;
  int exitcode, i, m, n, ngj;

  f = fopen(fn, "rb");
  if (!f)
    return -1;
  if (fread(&w->np, 2 * sizeof(int) + sizeof(float), 1, f) != 1) {
    exitcode = -2;
    goto err1;
  }
  w->n = (int *) calloc(3 * w->nfun, sizeof(int));
  if (!w->n) {
    exitcode = -3;
    goto err1;
  }
  w->nj = w->n + w->nfun;
  w->i = w->nj + w->nfun;
  if (fread(w->n, sizeof(int), w->nfun, f) != w->nfun) {
    exitcode = -4;
    goto err2;
  }
  ngj = (nge + 1) / reduc - 1;
  m = 0;
  for (i = 0; i < w->nfun; i++) {
    w->i[i] = m;
    m += w->n[i];
    n = (w->n[i] + 1) / reduc - 1;
    w->nj[i] = n > ngj ? ngj : n;
    if (w->n[i] > nge) w->n[i] = nge;
  }
  w->s = (float *) calloc(m, sizeof(float));
  if (!w->s) {
    exitcode = -5;
    goto err2;
  }
  if (fread(w->s, sizeof(float), m, f) != m) {
    exitcode = -6;
    goto err3;
  }
  fclose(f);
  return 0;
 err3:
  free(w->s);
 err2:
  free(w->n);
 err1:
  fclose(f);
  return exitcode;
}

void avgw_del(avgw_t *w)
{
  free(w->n);
  free(w->s);
}

#define max(a, b) ((a) > (b) ? a : b)
#define min(a, b) ((a) < (b) ? a : b)

void avgw_mul_dbl(const avgw_t *w, int ngrid, int natv, int natu, const double *x, double *r)
{
  int np, s, i0, n, u, v, j, i, l, k, incu;
  double a;

  incu = ngrid * natv;
  np = incu * natu;


  cblas_dcopy(np, x, 1, r, 1);
  for (s = 0; s < w->nfun; s++) {
    i0 = w->i[s];
    n = w->n[s];

    u = (int) (sqrt(2.0 * s + 0.25) + 0.5);
    j = s - u * (u - 1) / 2;

    l = u * incu;
    k = j * incu;

    for (v = 0; v < natv; v++) {
      for (i = 0; i < n; i++) {
        a = (double) w->s[i0+i];
        r[l+i] += a * x[k+i];
        r[k+i] += a * x[l+i];
      }
      l += ngrid;
      k += ngrid;
    }
  }
}

void avgw_mul_flt(const avgw_t *w, int incw, int ngrid, int natv, int natu, const float *x, float *r)
{
  int np, s, i0, n, u, v, j, i, l, k, t, incu;
  float a;

  incu = ngrid * natv;
  np = incu * natu;

  cblas_scopy(np, x, 1, r, 1);
  for (s = 0; s < w->nfun; s++) {
    i0 = w->i[s];
    n = w->nj[s];

    u = (int) (sqrt(2.0 * s + 0.25) + 0.5);
    j = s - u * (u - 1) / 2;

    l = u * incu;
    k = j * incu;

    for (v = 0; v < natv; v++) {
      t = incw - 1;
      for (i = 0; i < n; i++) {
        a = w->s[i0 + t];
        r[l+i] += a * x[k+i];
        r[k+i] += a * x[l+i];
        t += incw;
      }
      l += ngrid;
      k += ngrid;
    }
  }
}
