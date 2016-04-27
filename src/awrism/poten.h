#ifndef __RISM_POTEN_H
#define __RISM_POTEN_H

#include "grid.h"
#include "water.h"
#include "mol.h"

struct POTENUV
{
  int ngrid;
  int *atyp;
  double *u;
  double *asympr;
  double *asympk;
  double k;
};
typedef struct POTENUV potenuv_t;


int poten_mk(const grid_t *g, const water_t *w, const mol_t *m, double ngalpha, potenuv_t *e);
void poten_del(potenuv_t *e);

void uljuv(const grid_t *g, const water_t *w, const mol_t *m, double *u);
void ucoluv(const grid_t *g, const water_t *w, const mol_t *m, double *u);

#endif //__RISM_POTEN_H
