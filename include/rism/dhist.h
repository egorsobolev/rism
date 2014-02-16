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
  int f0;
  int fn;
  int hst0;
  int hstn;
  int *lm;
  int *ld;
  unsigned *hst;
};

typedef struct DHIST dhist_t;

int dhist_init(int np, double dr, const lavg_t *s, int n, int p, dhist_t *d);
void dhist_free(dhist_t *d);
void dhist_setlocal(int nfun, int n, int p, dhist_t *d);
void dhist_update(const float *l, dhist_t *d);

#ifdef MPI
#include <mpi.h>

int dhist_process_l(dhist_t *d, int nfrm, MPI_File in);
int dhist_read_hdr(dhist_t *d, int n, MPI_File f);
int dhist_write_hdr(dhist_t *d, MPI_File f);
int dhist_write_hist(dhist_t *d, MPI_File f);
#else
int dhist_process_l(dhist_t *d, int nfrm, FILE *in);
int dhist_read_hdr(dhist_t *d, int n, FILE *f);
int dhist_write_hdr(dhist_t *d, FILE *f);
int dhist_write_hist(dhist_t *d, FILE *f);
#endif

#endif /* __RISM_DHIST_H */
