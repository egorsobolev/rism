#ifndef __RISM_TDYN_H
#define __RISM_TDYN_H

#include "eqoz.h"

double etot(rismaw_t *r, const double *huv, const double *cuv);
double musc_hnc(rismaw_t *r, const double *huv, const double *cuv);
double musc_plhnc(rismaw_t *r, const double *huv, const double *cuv);
double mugf(rismaw_t *r, const double *huv, const double *cuv);

#endif /* __RISM_TDYN_H */
