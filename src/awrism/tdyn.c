#include "tdyn.h"

#include <math.h>
#include <stdio.h>

double etot(rismaw_t *r, const double *huv, const double *cuv)
{
	int i, j, u, v, l, incu;
	double a;
	incu = r->natv * r->ge.n;
	a = 0.0;
	l = 0;
	for (u = 0; u < r->natu; u++) {
		j = r->puv.atyp[u] * incu;
		for (v = 0; v < r->natv; v++) {
			for (i = 1; i <= r->ge.n; i++) {
				a += (huv[l] + 1.0) * r->puv.u[j + i - 1] * i * i * r->v.n[v];
				l++;
			}
			j += r->puv.ngrid;
		}
	}
	return r->puv.k * a;
}

double musc_hnc(rismaw_t *r, const double *huv, const double *cuv)
{
	int i, u, v, l, incu;
	double a;
	incu = r->natv * r->ge.n;
	a = 0.0;
	l = 0;
	for (u = 0; u < r->natu; u++) {
		for (v = 0; v < r->natv; v++) {
			for (i = 1; i <= r->ge.n; i++) {
				a += (0.5 * huv[l] * (huv[l] - cuv[l]) - cuv[l]) * i * i * r->v.n[v];
				l++;
			}
		}
	}
	return r->puv.k * a;
}

double musc_plhnc(rismaw_t *r, const double *huv, const double *cuv)
{
	int i, u, v, l, incu;
	double a;
	incu = r->natv * r->ge.n;
	a = 0.0;
	l = 0;
	for (u = 0; u < r->natu; u++) {
		for (v = 0; v < r->natv; v++) {
			for (i = 1; i <= r->ge.n; i++) {
				a += (0.5 * huv[l] * ((huv[l] < 0.0) * huv[l] - cuv[l]) - cuv[l]) * i * i * r->v.n[v];
				l++;
			}
		}
	}
	return r->puv.k * a;
}

double mugf(rismaw_t *r, const double *huv, const double *cuv)
{
	int i, u, v, l, incu;
	double a;
	incu = r->natv * r->ge.n;
	a = 0.0;
	l = 0;
	for (u = 0; u < r->natu; u++) {
		for (v = 0; v < r->natv; v++) {
			for (i = 1; i <= r->ge.n; i++) {
				a -= (0.5 * huv[l] * cuv[l] + cuv[l]) * i * i * r->v.n[v];
				l++;
			}
		}
	}
	return r->puv.k * a;
}
