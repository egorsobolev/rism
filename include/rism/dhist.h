#ifndef __RISM_DHIST_H
#define __RISM_DHIST_H

#include <rism/lavg.h>

struct DHIST
{
  double dr;
  int np;
  int n;
  int nfun;
  int nfrm;
  int *lm;
  int *ld;
  unsigned *hst;
};

typedef struct DHIST dhist_t;

int dhist_init(int np, double dr, const lavg_t *s, dhist_t *d);
void dhist_free(dhist_t *d);
void dhist_update(const float *l, dhist_t *d);
int dhist_write_hdr(const dhist_t *d, FILE *f);
int dhist_read_hdr(dhist_t *d, FILE *f);
int dhist_write_hist(const dhist_t *d, FILE *f);

#endif /* __RISM_DHIST_H */
