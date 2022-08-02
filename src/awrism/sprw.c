#include "sprw.h"

#include <stdlib.h>
#include <math.h>
#include <cblas.h>

int sprw_mk(mol_t *m, double prec, double r, int ngl, int ngs, double dk, sprw_t *w)
{
  double *l;
  int i, j, k, nd, nl, ns, a, b, c, n;
  double dx, dy, dz, d, t;

  w->nfun = m->natom * (m->natom - 1) / 2;

  if (m->natom < 2)
    return 0;

  l = (double *) malloc(w->nfun * sizeof(double));
  if (!l) {
    return -1;
  }
  w->i = (int *) malloc(w->nfun * sizeof(double));
  if (!w->i) {
    free(l);
    return -2;
  }

  dx = m->x[0] - m->x[1];
  dy = m->y[0] - m->y[1];
  dz = m->z[0] - m->z[1];
  l[0] = round(sqrt(dx * dx + dy * dy + dz * dz) / prec) * prec;
  w->nl = (l[0] < 3.0);
  nd = 1;
  w->i[0] = 0;

  k = 0;
  for (i = 2; i < m->natom; i++)
    for (j = 0; j < i; j++) {
      dx = m->x[i] - m->x[j];
      dy = m->y[i] - m->y[j];
      dz = m->z[i] - m->z[j];
      d = round(sqrt(dx * dx + dy * dy + dz * dz) / prec) * prec;
      
      a = 0;
      b = nd;
      while ((b - a) > 1) {
	c = a + b / 2;
	if (d > l[c]) a = c;
	else if (d < l[c]) b = c;
	else break;
      }

      if (l[c] == d) {
	w->i[k] = c;
      } else {
	if (d < 3.0)
	  w->nl++;
	if (c < nd && l[c] < d)
	  c++;
	w->i[k] = c;
	while (c < nd) {
	  t = l[c];
	  l[c] = d;
	  d = t;
	  c++;
	}
	l[nd] = d;
	nd++;
      }
      k++;
    }

  ns = nd - w->nl;
  w->s = (float *) malloc((ns * ngs + w->nl * ngl) * sizeof(float));
  if (!w->s) {
    free(l);
    free(w->i);
    return -3;
  }

  k = 0;
  for (j = 0; j < nd; j++) {
    d = l[j] * dk;
    n = j < nl ? ngl : ngs;
    for (i = 0; i < n; i++) {
      w->s[k] = (float) (sin(i * d) / (i * d));
      k++;
    }
  }
  w->p = (ngl - ngs) * w->nl;
  w->ngl = ngl;
  w->ngs = ngs;
  
  free(l);

  return 0;
}

void sprw_del(sprw_t *w)
{
  free(w->s);
  free(w->i);
}

#define max(a, b) ((a) > (b) ? a : b)
#define min(a, b) ((a) < (b) ? a : b)

void sprw_mul_dbl(const sprw_t *w, int ngrid, int natv, int natu, const double *x, double *r)
{
  int np, s, i0, n, u, v, j, i, l, k, incu;
  double a;

  incu = ngrid * natv;
  np = incu * natu;

  cblas_dcopy(np, x, 1, r, 1);
  for (s = 0; s < w->nfun; s++) {
    i0 = w->i[s];
    if (i0 < w->nl) {
      n = w->ngl;
      i0 *= w->ngl;
    } else {
      n = w->ngs;
      i0 = w->ngs * i0 + (w->ngl - w->ngs) * w->nl;
    }
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

void sprw_mul_flt(const sprw_t *w, int incw, int ngrid, int natv, int natu, const float *x, float *r)
{
  int np, s, i0, n, u, v, j, i, l, k, t, incu;
  float a;

  incu = ngrid * natv;
  np = incu * natu;

  cblas_scopy(np, x, 1, r, 1);
  for (s = 0; s < w->nfun; s++) {
    i0 = w->i[s];
    if (i0 < w->nl) {
      n = w->ngl;
      i0 *= w->ngl;
    } else {
      n = w->ngs;
      i0 = w->ngs * i0 + (w->ngl - w->ngs) * w->nl;
    }
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

