#include "closure.h"
#include "eqoz.h"

#include <math.h>

void hnc(rismaw_t *p, const double *tuv, double *cuv, float *dcdt, double *en)
{
  int i, j, l, k, u, v, incu;

  incu = p->natv * p->ge.n;

  l = 0;
  k = 0;
  *en = 0;
  for (u = 0; u < p->natu; u++) {
    j = p->puv.atyp[u] * incu;
    for (v = 0; v < p->natv; v++) {
      for (i = 0; i < p->gj.n; i++) {
	cuv[l] = exp(tuv[l] - p->puv.u[j + i]) - 1.0;
	dcdt[k] = (float) cuv[l];
	*en += (cuv[l] + 1.0) * p->puv.u[j + i] * (i + 1) * (i + 1) * p->v.n[v]; 
	cuv[l] -= tuv[l];
	l++;
	k++;
      }
      for (; i < p->ge.n; i++) {
	cuv[l] = exp(tuv[l] - p->puv.u[j + i]) - 1.0;
	*en += (cuv[l] + 1.0) * p->puv.u[j + i] * (i + 1) * (i + 1) * p->v.n[v]; 
	cuv[l] -= tuv[l];
	l++;
      }
      j += p->puv.ngrid;
    }
  }
  *en *= p->puv.k;
}
void hnc_c(rismaw_t *p, const double *tuv, double *cuv)
{
  int i, j, l, u, incu;

  incu = p->natv * p->ge.n;

  l = 0;
  for (u = 0; u < p->natu; u++) {
    j = p->puv.atyp[u] * incu;
    for (i = 0; i < incu; i++) {
      cuv[l] = exp(tuv[l] - p->puv.u[j + i]) - tuv[l] - 1.0;
      l++;
    }
  }
}

void plhnc(rismaw_t *p, const double *tuv, double *cuv, float *dcdt, double *en)
{
  int i, j, l, k, u, v, incu;
  double a;

  incu = p->natv * p->ge.n;

  l = 0;
  k = 0;
  *en = 0;
  for (u = 0; u < p->natu; u++) {
    j = p->puv.atyp[u] * incu;
    for (v = 0; v < p->natv; v++) {
      for (i = 0; i < p->gj.n; i++) {
	a = tuv[l] - p->puv.u[j + i];
	cuv[l] = exp((a < 0.0) * a) - 1.0;
	dcdt[k] = (float) cuv[l];
	cuv[l] += (a >= 0.0) * a;
	*en += (cuv[l] + 1.0) * p->puv.u[j + i] * (i + 1) * (i + 1) * p->v.n[v]; 
	cuv[l] -= tuv[l];
	l++;
	k++;
      }
      for (; i < p->ge.n; i++) {
	a = tuv[l] - p->puv.u[j + i];
	cuv[l] = exp((a < 0.0) * a) - 1.0 + (a >= 0.0) * a;
	*en += (cuv[l] + 1.0) * p->puv.u[j + i] * (i + 1) * (i + 1) * p->v.n[v];	
	cuv[l] -= tuv[l];
	l++;
      }
      j += p->puv.ngrid;
    }
  }
  *en *= p->puv.k;
}

void plhnc_c(rismaw_t *p, const double *tuv, double *cuv)
{
  int i, j, l, u, incu;
  double a;

  incu = p->natv * p->ge.n;

  l = 0;
  for (u = 0; u < p->natu; u++) {
    j = p->puv.atyp[u] * incu;
    for (i = 0; i < incu; i++) {
      a = tuv[l] - p->puv.u[j + i];
      cuv[l] = exp((a < 0.0) * a) - 1.0 + (a >= 0.0) * a - tuv[l];
      l++;
    }
  }
}
