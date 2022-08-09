#include "asymp.h"

#include <math.h>

void asymp(const grid_t *g, const water_t *w, const mol_t *m, double alpha, double *asympr, double *asympk)
{
  int u, v, i, j, np;
  double er, ex, quv, a, b, thi, r, k;
  a = -18.2223 * 18.2223 / w->m.t;
  a = -332.0 / w->m.t;
  b = 4.0 * M_PI;
  thi = -0.25 / (alpha * alpha);
  np = g->n - 1;

  for (i = 1; i < g->n; i++) {
    j = i - 1;
    r = i * g->dr;
    k = i * g->dk;

    er = erf(alpha * r) / r;
    ex = b * exp(thi * k * k) / k;

    for (u = 0; u < m->ntype; u++)
      for (v = 0; v < w->natom; v++) {
        quv = a * w->q[v] * m->ff.q[u];
        asympr[j] = quv * er;
        asympk[j] = quv * ex;
        j += np;
      }
  }
}
