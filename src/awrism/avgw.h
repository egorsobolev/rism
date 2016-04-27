#ifndef __RISM_AVGW_H
#define __RISM_AVGW_H

struct AVGW
{
  int np;
  float dr;
  int nfun;
  int *n;
  int *i;
  int *nj;
  float *s;
};

typedef struct AVGW avgw_t;

int avgw_read(const char *fn, avgw_t *w, int nge, int reduc);
void avgw_del(avgw_t *w);

void avgw_mul_dbl(const avgw_t *w, int ngrid, int natv, int natu, const double *x, double *r);
void avgw_mul_flt(const avgw_t *w, int incw, int ngrid, int natv, int natu, const float *x, float *r);

#endif //__RISM_AVGW_H
