#ifndef __RISM_LAVG_H
#define __RISM_LAVG_H

#include <stdio.h>

struct LAVG {
  int natm;
  int nfun;
  int nfrm;
  double *mx;
  double *mn;
  double *mean;
  double *std;
};

typedef struct LAVG lavg_t;

int lavg_init(int natm, lavg_t *s);
void lavg_free(lavg_t *s);
void lavg_update(const float *x, float *l, lavg_t *s);
int lavg_finish(lavg_t *s);

#ifdef MPI
#include <mpi.h>

int lavg_readhdr(lavg_t *s, MPI_File f);
int lavg_writehdr(const lavg_t *s, MPI_File f);
#else
int lavg_readhdr(lavg_t *s, FILE *f);
int lavg_writehdr(const lavg_t *s, FILE *f);
#endif

#endif /* __RISM_LAVG_H */
