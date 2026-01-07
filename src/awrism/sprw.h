#ifndef __RISM_SPRW_H
#define __RISM_SPRW_H

#include "mol.h"

struct SPRW
{
	int nfun;
	/*
	int *n;
	int *i0;
	int *nj;
	*/
	int ngs;
	int ngl;
	int nl;
	int p;
	int *i;
	float *s;
};
typedef struct SPRW sprw_t;

int sprw_mk(mol_t *m, double prec, double r, int ngl, int ngs, double dk, sprw_t *w);
void sprw_del(sprw_t *w);

void sprw_mul_dbl(const sprw_t *w, int ngrid, int natv, int natu, const double *x, double *r);
void sprw_mul_flt(const sprw_t *w, int incw, int ngrid, int natv, int natu, const float *x, float *r);

#endif /* __RISM_SPRW_H */
