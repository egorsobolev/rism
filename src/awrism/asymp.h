#ifndef __RISM_ASYMP_H
#define __RISM_ASYMP_H

#include "grid.h"
#include "water.h"
#include "mol.h"

void asymp(const grid_t *g, const water_t *w, const mol_t *m, double alpha, double *asympr, double *asympk);

#endif /* __RISM_ASYMP_H */
