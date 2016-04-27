#ifndef __RISM_EQOZ_H
#define __RISM_EQOZ_H

#include "avgw.h"
#include "poten.h"

struct SOLV
{
  int ngrid;
  double *xvv;
  double *symc;
  int *n;
  double t;
  double rho;
};
typedef struct SOLV solv_t;

#include "grid.h"

typedef struct RISM_IDL_AW rismaw_t;
typedef void closure(rismaw_t *, const double *, double *, float *, double *);
typedef void closure_c(rismaw_t *, const double *, double *);

struct RISM_IDL_PARM
{
  int reduc;
  double ngalpha;
};
typedef struct RISM_IDL_PARM rism_parm_t;

struct RISM_IDL_AW
{
  int natu;
  int natv;
  int nfun;

  solv_t v;
  grid_eq_t ge;
  avgw_t wuu;
  potenuv_t puv;
  closure *closure;
  closure_c *closure_c;

  grid_jac_t gj;
  float *dcdt;

  int reduc;
  double ngalpha;
};


int rismaw_eq(void *p, const double *tuv, double *d, double *en);
int rismaw_Jx(void *p, const float *x, float *r);
void rismaw_getx(const void *p, const double *d, float *x);
void rismaw_putx(const void *p, const float *x, double *d);

#endif //__RISM_EQOZ_H
