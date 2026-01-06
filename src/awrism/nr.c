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
#define flush	fflush(stdin)

#define sarm    (0.618033988749895)


/* + 5 * n * sizeof(float)*/
int bicgstab(int N, const float *b, float *x, Jx_func *Jx, void *eq_data, float *tol, int *it);

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
  double *d, *s, *tN, norms, eta, lambda, err, errN, v, en1, en2, dE, errv;
  float *b, *x, eps;
  int n, m, i, j, flag, nlit, maxlit, l, narm;
  double lntm, nrtm, lntmi, nrtm_cpu;
  fterm_ew_t fterm;

  nrtm_cpu = clock() / (double) CLOCKS_PER_SEC;
  nrtm = walltime() * 1e-6;
  lntm = .0;

  ftew_setdef(&fterm, *rtol);
  maxlit = eq->nrprm.lmaxi;
  narm = eq->nrprm.narm;
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
  v = sqrt(n);

  b = (float *) malloc(2 * m * sizeof(float) + 3 * n * sizeof(double));
  if (!b) {
    return -1;
  }
  x = b + m;
  d = (double *) (x + m);
  s = d + n;
  tN = s + n;

  eq->Z(eq->p, t, d, &en1);
  err = cblas_dnrm2(n, d, 1) / v;

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

    eq->getb(eq->p, d, b);
    memset(x, 0, m * sizeof(float));

    lntmi = walltime() * 1e-6;
    flag = bicgstab(m, b, x, eq->Jx, eq->p, &eps, &nlit);
    lntmi = walltime() * 1e-6 - lntmi;
    lntm += lntmi;

    if (flag) {
      /* no convergence of linear solver */
      return 1;
    }
    norms = (double) eps;

    cblas_dcopy(n, d, 1, s, 1);
    cblas_sscal(m, -1.0, x, 1);
    eq->putx(eq->p, x, s);

    cblas_dcopy(n, t, 1, tN, 1);
    cblas_daxpy(n, 1.0, s, 1, tN, 1);
    en2 = en1;
    eq->Z(eq->p, tN, d, &en1);
    errN = cblas_dnrm2(n, d, 1) / v;

    j = 0;
    lambda = 1.0;
    while (j < narm && errN > err) {
      lambda *= sarm;
      cblas_dcopy(n, t, 1, tN, 1);
      cblas_daxpy(n, lambda, s, 1, tN, 1);
      eq->Z(eq->p, tN, d, &en1);
      errN = cblas_dnrm2(n, d, 1) / v;
      ++j;
    }

    if (j >= narm) {
      /* no convergence of armigo*/
      return 2;
    }
    cblas_dcopy(n, tN, 1, t, 1);
    err = errN;
    dE = en2 - en1;
    errv = lambda * cblas_dnrm2(n, s, 1) / v;

    print("%5d %7d %1d %8.1f %1d %5.2f %10.6f %8.1e %8.1e %8.1e\n", i+1, nlit, flag, lntmi, j, eta, norms, err, dE, errv);
    flush;
    ++i;
  }

  free(b);

  *maxit = i;
  *rtol = err;
  *etol = abs(dE);

  nrtm_cpu = clock() / (double) CLOCKS_PER_SEC - nrtm_cpu;
  nrtm = walltime() * 1e-6 - nrtm;
  print("NR: %.1lfs, BiCGStab: %.1lfs, CPU usage: %.1lf\n", nrtm, lntm, nrtm_cpu / nrtm);
  print("RMSE(Z[%d]) = %.1e\n", i, err);

  return 0; /* OK */
}
