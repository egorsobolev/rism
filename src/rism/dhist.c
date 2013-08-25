#include <rism/dhist.h>
#include <rism/lavg.h>

#include <stdlib.h>
#include <memory.h>
#include <stdio.h>

int dhist_init(int np, double dr, const lavg_t *s, dhist_t *d)
{
  int m, i;
  double a;
  d->np = np;
  d->dr = dr;
  d->nfun = s->nfun;
  d->nfrm = s->nfrm;

  m = 2 * d->nfun;
  d->lm = (int *) calloc(m, sizeof(int));
  if (!d->lm)
    return -1;
  d->ld = d->lm + d->nfun;

  d->n = 0;
  for (i = 0; i < s->nfun; i++) {
    d->lm[i] = (int) (s->mn[i] / dr + .5);
    d->ld[i] = (int) (s->mx[i] / dr + .5) - d->lm[i] + 1;
    d->n += d->ld[i];
  }
  d->hst = (unsigned *) calloc(d->n, sizeof(unsigned));
  if (!d->hst) {
    free(d->lm);
    return -1;
  }

  return 0;
}

void dhist_free(dhist_t *d)
{
  free(d->lm);
  free(d->hst);
}

void dhist_update(const float *l, dhist_t *d)
{
  int k, i;
  k = 0;
  for (i = 0; i < d->nfun; i++) {
    d->hst[k + (int) (l[i] / d->dr + 0.5) - d->lm[i]]++;
    k += d->ld[i];
  }
}

int dhist_write_hdr(const dhist_t *d, FILE *f)
{
  int m;
  if (fwrite(&d->dr, sizeof(double), 1, f) != 1)
    return -1;
  if (fwrite(&d->np, sizeof(int), 4, f) != 4)
    return -1;
  m = 2 * d->nfun;
  if (fwrite(d->lm, sizeof(int), m, f) != m)
    return -1;
  return 0;
}

int dhist_read_hdr(dhist_t *d, FILE *f)
{
  int m;
  if (fread(&d->dr, sizeof(double), 1, f) != 1)
    return -1;
  if (fread(&d->np, sizeof(int), 4, f) != 4)
    return -1;

  m = 2 * d->nfun;
  d->lm = (int *) calloc(m, sizeof(int));
  if (!d->lm)
    return -1;
  d->ld = d->lm + d->nfun;
  if (fread(d->lm, sizeof(int), m, f) != m) {
    free(d->lm);
    return -1;
  }
  d->hst = (unsigned *) calloc(d->np, sizeof(unsigned));
  if (!d->hst) {
    free(d->lm);
    return -1;
  }
  return 0;
}

int dhist_write_hist(const dhist_t *d, FILE *f) 
{
  return fwrite(d->hst, sizeof(unsigned), d->n, f) != d->n ? -1: 0;
}
