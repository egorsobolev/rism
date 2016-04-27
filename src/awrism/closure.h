#ifndef __RISM_CLOSURE_H
#define __RISM_CLOSURE_H

#include "eqoz.h"

void hnc(rismaw_t *p, const double *tuv, double *cuv, float *dcdt, double *en);
void hnc_c(rismaw_t *p, const double *tuv, double *cuv);
void plhnc(rismaw_t *p, const double *tuv, double *cuv, float *dcdt, double *en);
void plhnc_c(rismaw_t *p, const double *tuv, double *cuv);

#endif //__RISM_CLOSURE_H
