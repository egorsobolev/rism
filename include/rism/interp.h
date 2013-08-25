#ifndef __RISM_INTERP_H
#define __RISM_InTERP_H

void sspline(int n1, const float *x, const float *y, float yp1, float ypn, float *y2, float *w);
void ssplint(int n, const float *xa, const float *ya, const float *y2a, int m, float *xx, float *yy);
void sspline_uni(int n1, const float *y, float yp1, float ypn, float *y2, float *w);
void ssplint_uni(int n, const float *ya, const float *y2a, float m, float x0, float h, float *yy);

#endif /* __RISM_INTERP_H */
