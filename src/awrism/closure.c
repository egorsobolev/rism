#include "closure.h"
#include "eqoz.h"

#include <math.h>

void hnc(rismaw_t *p, const double *tuv, double *cuv, float *dcdt, double *en)
{
	int i, j, l, k, u, v, incu, t;
	double eni;

	incu = p->natv * p->ge.n;

	eni = 0;
	#pragma omp single nowait
	*en = 0;
	#pragma omp for collapse(2)
	for (u = 0; u < p->natu; u++) {
		for (v = 0; v < p->natv; v++) {
			j = p->puv.atyp[u] * incu + v * p->puv.ngrid;
			t = (v + p->natv * u);
			for (i = 0; i < p->gj.n; i++) {
				k = i + p->gj.n * t;
				l = i + p->ge.n * t;
				cuv[l] = exp(tuv[l] - p->puv.u[j + i]) - 1.0;
				dcdt[k] = (float) cuv[l];
				eni += (cuv[l] + 1.0) * p->puv.u[j + i] * (i + 1) * (i + 1) * p->v.n[v];
				cuv[l] -= tuv[l];
			}
			for (; i < p->ge.n; i++) {
				l = i + p->ge.n * t;
				cuv[l] = exp(tuv[l] - p->puv.u[j + i]) - 1.0;
				eni += (cuv[l] + 1.0) * p->puv.u[j + i] * (i + 1) * (i + 1) * p->v.n[v];
				cuv[l] -= tuv[l];
			}
		}
	}
	#pragma atomic
	*en += p->puv.k * eni;
}
void hnc_c(rismaw_t *p, const double *tuv, double *cuv)
{
	int i, j, l, u, incu;

	incu = p->natv * p->ge.n;

	#pragma omp for collapse(2)
	for (u = 0; u < p->natu; u++) {
		for (i = 0; i < incu; i++) {
			j = p->puv.atyp[u] * incu;
			l = i + u * incu;
			cuv[l] = exp(tuv[l] - p->puv.u[j + i]) - tuv[l] - 1.0;
		}
	}
}

void plhnc(rismaw_t *p, const double *tuv, double *cuv, float *dcdt, double *en)
{
	int i, j, l, k, u, v, incu, t;
	double a, eni;

	incu = p->natv * p->ge.n;

	eni = 0;
	#pragma omp single nowait
	*en = 0;
	#pragma omp for collapse(2)
	for (u = 0; u < p->natu; u++) {
		for (v = 0; v < p->natv; v++) {
			j = p->puv.atyp[u] * incu + v * p->puv.ngrid;
			t = (v + p->natv * u);
			for (i = 0; i < p->gj.n; i++) {
				k = i + p->gj.n * t;
				l = i + p->ge.n * t;
				a = tuv[l] - p->puv.u[j + i];
				cuv[l] = exp((a < 0.0) * a) - 1.0;
				dcdt[k] = (float) cuv[l];
				cuv[l] += (a >= 0.0) * a;
				eni += (cuv[l] + 1.0) * p->puv.u[j + i] * (i + 1) * (i + 1) * p->v.n[v];
				cuv[l] -= tuv[l];
			}
			for (; i < p->ge.n; i++) {
				l = i + p->ge.n * t;
				a = tuv[l] - p->puv.u[j + i];
				cuv[l] = exp((a < 0.0) * a) - 1.0 + (a >= 0.0) * a;
				eni += (cuv[l] + 1.0) * p->puv.u[j + i] * (i + 1) * (i + 1) * p->v.n[v];
				cuv[l] -= tuv[l];
			}
		}
	}
	#pragma omp atomic
	*en += p->puv.k * eni;
}

void plhnc_c(rismaw_t *p, const double *tuv, double *cuv)
{
	int i, j, l, u, incu;
	double a;

	incu = p->natv * p->ge.n;

	#pragma omp for collapse(2)
	for (u = 0; u < p->natu; u++) {
		for (i = 0; i < incu; i++) {
			j = p->puv.atyp[u] * incu;
			l = i + u * incu;
			a = tuv[l] - p->puv.u[j + i];
			cuv[l] = exp((a < 0.0) * a) - 1.0 + (a >= 0.0) * a - tuv[l];
		}
	}
}
