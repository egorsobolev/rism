#include "eqoz.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <fft.h>

void rismaw_getx(const void *prm, const double *d, float *x)
{
	int i, j;
	const rismaw_t *p = (rismaw_t *) prm;

	#pragma omp for collapse(2)
	for (i = 0; i < p->nfun; i++)
		for (j = 0; j < p->gj.n; j++)
			x[i * p->gj.n + j] = (float) d[i * p->ge.n + j];
}

void rismaw_putx(const void *prm, const float *x, double *d)
{
	int i, j;
	const rismaw_t *p = (rismaw_t *) prm;

	#pragma omp for collapse(2)
	for (i = 0; i < p->nfun; i++)
		for (j = 0; j < p->gj.n; j++)
			d[i * p->ge.n + j] = (double) x[i * p->gj.n + j];
}

#define max(a, b) ((a) > (b) ? a : b)
#define min(a, b) ((a) < (b) ? a : b)

int rismaw_eq(void *prm, const double *tuv, double *d, double *en)
{
	double *r;
	int i, j, u, v, k, l, s, np, incu, i0, n, t, vng1;
	double a;

	rismaw_t *p = (rismaw_t *) prm;
	grid_t *g = &p->ge;

	incu = p->natv * g->n;
	np = p->natu * incu;
	vng1 = p->v.ngrid + 1;

	r = p->Z_data;

	/* cuv = C[tuv(r)] */
	p->closure(p, tuv, r, p->dcdt, en);

	/* r = fft(cuv) */
	n = min(p->puv.ngrid, g->n);
	#pragma omp for collapse(2)
	for (u = 0; u < p->natu; u++) {
		for (v = 0; v < p->natv; v++) {
			l = (v + u * p->natv) * g->n;
			k = p->puv.atyp[u] * p->natv * p->puv.ngrid + v * p->puv.ngrid;
			for (i = 0; i < n; i++)
				r[l + i] = (r[l + i] - p->puv.asympr[k + i]) * (i + 1);
			for (; i < g->n; i++)
				r[l + i] = 0.0;

			fft_dst(1, &g->n, 1, r + l, r + l, 1.0, 0);

			for (i = 0; i < n; i++)
				r[l + i] = (g->f * r[l + i] + p->puv.asympk[k + i]) * p->v.symc[v];
			for (; i < g->n; i++)
				r[l + i] = 0.0;
		}
	}

	/* d = wuu * r = wuu * fft(cuv) */
	avgw_mul_dbl(&p->wuu, g->n, p->natv, p->natu, r, d);

	/* r = d * xvv - r = wuu * fft(cuv) * xvv - fft(cuv) */
	#pragma omp for
	for (i = 0; i < np; i++)
		r[i] = -r[i];
	n = min(p->v.ngrid - 1, g->n);
	for (v = 0; v < p->natv; v++) {
		/* side elements */
		i0 = v * (v + 1) / 2 * vng1 + 1;
		for (j = 0; j < v; j++) {
			l = v * g->n;
			k = j * g->n;
			#pragma omp for collapse(2)
			for (u = 0; u < p->natu; u++) {
				for (i = 0; i < n; i++) {
					t = u * incu + i;
					a = p->v.xvv[i0 + i + j * vng1];
					r[l + t] += a * d[k + t];
					r[k + t] += a * d[l + t];
				}
			}
		}
		/* diagonal */
		l = v * g->n;
		i0 += v * vng1;
		#pragma omp for collapse(2)
		for (u = 0; u < p->natu; u++) {
			for (i = 0; i < n; i++) {
				t = l + u * incu + i;
				r[t] += d[t] * p->v.xvv[i0 + i];
			}
		}
	}

	/* d = ifft(r) - tuv = ifft(wuu * fft(cuv) * xvv - fft(cuv)) - tuv */
	n = min(p->puv.ngrid, g->n);
	#pragma omp for collapse(2)
	for (u = 0; u < p->natu; u++) {
		for (v = 0; v < p->natv; v++) {
			l = (v + u * p->natv) * g->n;
			k = p->puv.atyp[u] * p->natv * p->puv.ngrid + v * p->puv.ngrid;
			for (i = 0; i < n; i++)
				d[l + i] = r[l + i] / p->v.symc[v] + p->puv.asympk[k + i];
			for (; i < g->n; i++)
				d[l + i] = 0.0;

			fft_dst(1, &g->n, 1, d + l, d + l, 1.0, 0);

			for (i = 0; i < n; i++)
				d[l + i] = g->b * d[l + i] / (i + 1) - p->puv.asympr[k + i] - tuv[l + i];

			for (; i < g->n; i++)
				d[l + i] = -tuv[l + i];

		}
	}

	return 0;
}

int rismaw_Jx(void *prm, const float *x, float *r)
{
	float *d;
	int i, j, u, v, k, l, s, np, incu, i0, n, t, vng1;
	float a;

	rismaw_t *p = (rismaw_t *) prm;
	grid_t *g = &p->gj;

	incu = p->natv * g->n;
	np = p->natu * incu;
	vng1 = p->v.ngrid + 1;

	d = p->Jx_data;

	/* Jx = fft(w * fft(dcdt * x) * xvv - fft(dcdt * x)) - x */

	/* r1 = fft(dcdt * x) */
	#pragma omp for collapse(2)
	for (u = 0; u < p->natu; u++) {
		for (v = 0; v < p->natv; v++) {
			l = (v + u * p->natv) * g->n;
			for (i = 0; i < g->n; i++)
				r[l + i] = p->dcdt[l + i] * x[l + i] * (i + 1);

			a = g->f * (float) p->v.symc[v];
			fftf_dst(1, &g->n, 1, r + l, r + l, a, 0);
		}
	}

	/* d = wuu * r = wuu * fft(cuv) */
	avgw_mul_flt(&p->wuu, p->reduc, g->n, p->natv, p->natu, r, d);

	/* r = d * xvv - r = wuu * fft(dcdt * x) * xvv - fft(dcdt * x) */
	#pragma omp for
	for (i = 0; i < np; i++)
		r[i] = -r[i];
	n = min(p->v.ngrid / p->reduc - 1, g->n);
	for (v = 0; v < p->natv; v++) {
		i0 = v * (v + 1) / 2 * (p->v.ngrid + 1);
		/* side elements */
		for (j = 0; j < v; j++) {
			l = v * g->n;
			k = j * g->n;
			#pragma omp for collapse(2)
			for (u = 0; u < p->natu; u++) {
				for (i = 0; i < n; i++) {
					t = u * incu + i;
					a = (float) p->v.xvv[i0 + (i + 1) * p->reduc + j * vng1];
					r[l + t] += a * d[k + t];
					r[k + t] += a * d[l + t];
				}
			}
		}
		/* diagonal */
		l = v * g->n;
		i0 += v * vng1;
		#pragma omp for collapse(2)
		for (u = 0; u < p->natu; u++) {
			for (i = 0; i < n; i++) {
				t = l + u * incu + i;
				r[t] += d[t] * (float) p->v.xvv[i0 + (i + 1) * p->reduc];
			}
		}
	}

	/* r = ifft(r3) - x = ifft(wuu * fft(dcdt * x) * xvv - fft(dcdt * x)) - x */
	#pragma omp for collapse(2)
	for (u = 0; u < p->natu; u++) {
		for (v = 0; v < p->natv; v++) {
			l = (v + u * p->natv) * g->n;
			for (i = 0; i < g->n; i++)
				r[l + i] /= p->v.symc[v];

			fftf_dst(1, &g->n, 1, r + l, r + l, 1.0, 0);

			for (i = 0; i < g->n; i++)
				r[l + i] = g->b * r[l + i] / (i + 1) - x[l + i];
		}
	}

	return 0;
}
