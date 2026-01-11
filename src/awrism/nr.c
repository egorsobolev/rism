#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include <cblas.h>
#include <time.h>

#include <rism/utils.h>
#include "nr.h"

#define print   printf
#define flush   fflush(stdin)

#define SCALE_ARM    (0.618033988749895)


/* + 5 * n * sizeof(float)*/
int bicgstab(int N, const float *b, float *x, Jx_func *Jx, void *eq_data, float *tol, int *it, float *workspace);

/* Eisenstat, Walker */
typedef struct
{
	double etaM;
	double eta;
	double gamma;
	double prec;
	double normdT;
} fterm_ew_t;

void ftew_setdef(fterm_ew_t *f, double tol)
{
	f->etaM = 0.9;
	f->eta = sqrt(f->etaM / 0.9);
	f->gamma = 0.9;
	f->prec = tol;
	f->normdT = 1.0;
}

double ftew_next(fterm_ew_t *f, double normd)
{
	double etaT;
	etaT = f->gamma * f->eta * f->eta;
	f->eta = normd / f->normdT;
	f->eta = f->gamma * f->eta * f->eta;
	if (etaT > .1 && etaT > f->eta)
		f->eta = etaT;
	if (f->etaM < f->eta)
		f->eta = f->etaM;
	etaT = f->prec / normd;
	if (etaT > f->eta)
		f->eta = etaT;
	f->normdT = normd;

	return f->eta;
}

/* 9 * n * sizeof(float) */
int nr(eq_t *eq, double *t, double *rtol, double *etol, int *maxit)
{
	double *d, *s, *tN, norms, eta, lambda, err, errN, sqrt_n, en1, en2, dE, errv;
	float *b, *x, eps;
	int n, m, i, j, flag, nlit, maxlit, l, narm, maxarm, fail;
	double sumd2;
	double lntm, nrtm, lntmi, nrtm_cpu;
	fterm_ew_t fterm;

	nrtm_cpu = clock() / (double) CLOCKS_PER_SEC;
	nrtm = walltime() * 1e-6;
	lntm = .0;

	ftew_setdef(&fterm, *rtol);
	maxlit = eq->nrprm.lmaxi;
	maxarm = eq->nrprm.narm;
	/*
	ln.natv = eq->solv.natv;
	ln.symc = eq->solv.symc;
	ln.lxvva = eq->solv.lxvva;
	ln.xvva = eq->solv.xvva;
	*/
	/*
		m = ln.natv * ln.natv * ln.lxvva;
		ln.xvva = malloc(m * sizeof(float));
		for (j = 0; j < m; ++j)
		ln.xvva[j] = (float) eq->solv.xvva[j];
	*/
	/*
	ln.indga = eq->solv.indga;
	ln.grid = &eq->lngr;
	m = ln.grid->nr * ln.natv;
	ln.dcdg = (float *) malloc(m * sizeof(float));
	*/

	n = eq->nZ;
	m = eq->nJx;
	/*
	n = eq->solv.natv * eq->grid.nr;
	*/
	sqrt_n = sqrt(n);

	b = eq->nr_solver_data;
	x = b + m;
	d = (double *) (x + m);
	s = d + n;
	tN = s + n;

	fail = 0;
	sumd2 = 0.0;
	#pragma omp parallel
	{
		int k;

		eq->Z(eq->p, t, d, &en1);

		#pragma omp for reduction(+:sumd2)
		for (k = 0; k < n; k++)
			sumd2 += d[k] * d[k];

		#pragma omp single nowait
		err = sqrt(sumd2) / sqrt_n;

		#pragma omp for
		for (k = 0; k < m; k++)
			x[k] = .0f;

		eq->getb(eq->p, d, b);

	}
	dE = DBL_MAX;

	print("RMSE(Z[0]) = %.1e\n", err);
	print("    #     Nit F      t,s M   eta   ||Jx-Z||  RMSE(Z)       dE  RMSE(v)\n");
	flush;
	i = 0;
	while (i <= *maxit && (err > *rtol || fabs(dE) > *etol)) {

		nlit = maxlit;
		eta = ftew_next(&fterm, err);
		eps = (float) eta;
		eps = 0.2;


		lntmi = walltime() * 1e-6;
		flag = bicgstab(m, b, x, eq->Jx, eq->p, &eps, &nlit, eq->linear_solver_data);
		lntmi = walltime() * 1e-6 - lntmi;
		lntm += lntmi;

		if (flag) {
			/* no convergence of linear solver */
			return 1;
		}
		norms = (double) eps;
		en2 = en1;

		#pragma omp parallel
		{
			int k, j;
			double errN, lambda;

			#pragma omp for
			for (k = 0; k < n; k++)
				s[k] = -d[k];

			eq->putx(eq->p, x, s);

			#pragma omp for
			for (k = 0; k < n; k++)
				tN[k] = t[k] - s[k];

			eq->Z(eq->p, tN, d, &en1);

			#pragma single nowait
			sumd2 = 0.0;
			#pragma omp for reduction(+:sumd2)
			for (k = 0; k < n; k++)
				sumd2 += d[k] * d[k];
			errN = sqrt(sumd2) / sqrt_n;

			j = 0;
			lambda = 1.0;
			while (j < maxarm && errN > err) {
				lambda *= SCALE_ARM;
				#pragma omp for
				for (k = 0; k < n; k++)
					tN[k] = t[k] - lambda * s[k];

				eq->Z(eq->p, tN, d, &en1);

				#pragma omp single nowait
				sumd2 = 0.0;
				#pragma omp for reduction(+:sumd2)
				for (k = 0; k < n; k++)
					sumd2 += d[k] * d[k];
				errN = sqrt(sumd2) / sqrt_n;
				#pragma omp single
				printf("* %d %lf %lf\n", j, errN, err);
				++j;
			}
			if (j >= maxarm) {
				/* no convergence of armigo*/
				#pragma omp single nowait
				fail = 1;
			} else {
				#pragma omp barrier
				#pragma omp single
				sumd2 = 0.0;
				#pragma omp for reduction(+:sumd2)
				for (k = 0; k < n; k++) {
					t[k] = tN[k];
					sumd2 += s[k] * s[k];
				}

				#pragma omp for
				for (k = 0; k < m; k++)
					x[k] = .0f;

				eq->getb(eq->p, d, b);

				#pragma omp single nowait
				errv = lambda * sqrt(sumd2) / sqrt_n;
				#pragma omp single nowait
				err = errN;
				#pragma omp single nowait
				dE = en2 - en1;
				#pragma omp single nowait
				narm = j;
			}
		}
		if (fail) {
			/* no convergence of armigo*/
			return 2;
		}
		print("%5d %7d %1d %8.1f %1d %5.2f %10.6f %8.1e %8.1e %8.1e\n", i+1, nlit, flag, lntmi, narm, eta, norms, err, dE, errv);
		flush;
		++i;
	}

	*maxit = i;
	*rtol = err;
	*etol = abs(dE);

	nrtm_cpu = clock() / (double) CLOCKS_PER_SEC - nrtm_cpu;
	nrtm = walltime() * 1e-6 - nrtm;
	print("NR: %.1lfs, BiCGStab: %.1lfs, CPU usage: %.1lf\n", nrtm, lntm, nrtm_cpu / nrtm);
	print("RMSE(Z[%d]) = %.1e\n", i, err);

	return 0; /* OK */
}
