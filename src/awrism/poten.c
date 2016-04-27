#include "poten.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cblas.h>

int poten_mk(const grid_t *g, const water_t *w, const mol_t *m, double ngalpha, potenuv_t *puv)
{
  int np;
  /*
  double a, b, c;

  a = 4.0 * M_PI * w->m.rho * g->dr * g->dr * g->dr;
  nfun = m->ntype * w->natom;
  */
  puv->ngrid = g->n - 1;
  np = m->ntype * w->natom * puv->ngrid;
  puv->u = (double *) malloc(3 * np * sizeof(double));
  if (!puv->u)
    return -1;
  puv->asympr = puv->u + np;
  puv->asympk = puv->asympr + np;

  memset(puv->u, 0, np * sizeof(double));
  ucoluv(g, w, m, puv->u);
  uljuv(g, w, m, puv->u);
  cblas_dscal(np, 1.0 / w->m.t, puv->u, 1);

  asymp(g, w, m, ngalpha, puv->asympr, puv->asympk);

  puv->atyp = m->atyp;

  puv->k = 4.0 * M_PI * g->dr * g->dr * g->dr * w->m.t * w->m.rho;

  return 0;
}

void poten_del(potenuv_t *puv)
{
  free(puv->u);
}

void uljuv(const grid_t *g, const water_t *w, const mol_t *m, double *phi)
{
  int u, v, i, j;
  double s, e, r6, a;
  a = 4.0;
  j = 0;
  for (u = 0; u < m->ntype; u++)
    for (v = 0; v < w->natom; v++) {
      s = 0.5 * (w->s[v] + m->ff.s[u]);
      e = a * sqrt(w->e[v] * m->ff.e[u]);
      for (i = 1; i < g->n; i++) {
	r6 = s / (i * g->dr);
	r6 = r6 * r6;
	r6 = r6 * r6 * r6;
	phi[j] += e * r6 * (r6 - 1.0);
	j++;
      }
    }
}

void ucoluv(const grid_t *g, const water_t *w, const mol_t *m, double *phi)
{
  int u, v, i, j;
  double quv, a;
  a = 18.2223 * 18.2223;
  a = 332.0;
  j = 0;
  for (u = 0; u < m->ntype; u++)
    for (v = 0; v < w->natom; v++) {
      quv = a * w->q[v] * m->ff.q[u];
      for (i = 1; i < g->n; i++) {
	phi[j] += quv / (i * g->dr);
	j++;
      }
    }
}
