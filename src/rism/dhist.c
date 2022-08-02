#include <rism/dhist.h>
#include <rism/lavg.h>
#include <rism/error.h>

#include <stdlib.h>
#include <memory.h>
#include <stdio.h>

#ifdef MPI
#include <mpi/mpi.h>
#endif

int dhist_init(int np, double dr, const lavg_t *s, int n, int p, dhist_t *d)
{
  int m, i;
  double a;
  d->np = np;
  d->dr = dr;
  d->nfun = s->nfun;
  d->nfrm = s->nfrm;

  dhist_setlocal(s->nfun, n, p, d);

  m = 2 * d->fn;
  d->lm = (int *) calloc(m, sizeof(int));
  if (!d->lm)
    return -1;
  d->ld = d->lm + d->fn;

  d->hst0 = 0;
  for (i = 0; i < d->f0; i++)
    d->hst0 += (int) (s->mx[i] / dr + .5) - (int) (s->mn[i] / dr + .5) + 1;
  
  d->hstn = 0;
  for (i = 0; i < d->fn; i++) {
    d->lm[i] = (int) (s->mn[d->f0 + i] / dr + .5);
    d->ld[i] = (int) (s->mx[d->f0 + i] / dr + .5) - d->lm[i] + 1;
    d->hstn += d->ld[i];
  }
  d->hst = (unsigned *) calloc(d->hstn, sizeof(unsigned));
  if (!d->hst) {
    free(d->lm);
    return -1;
  }

  d->n = d->hst0 + d->hstn;
  for (i = d->f0 + d->fn; i < d->nfun; i++)
    d->n += (int) (s->mx[i] / dr + .5) - (int) (s->mn[i] / dr + .5) + 1;

  return 0;
}

void dhist_setlocal(int nfun, int n, int p, dhist_t *d)
{
  int r;
  d->fn = nfun / n;
  r = nfun - d->fn * n;
  if (p < r) {
    d->fn += 1;
    d->f0 = d->fn * p;
  } else {
    d->f0 = d->fn * p + r;
  }
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
  for (i = 0; i < d->fn; i++) {
/*    printf("l[%d] = %f, i = %d, lm = %d\n", i, l[i],(int) (l[i] / d->dr + 0.5) - d->lm[i],d->lm[i]);*/
    d->hst[k + (int) (l[i] / d->dr + 0.5) - d->lm[i]]++;
    k += d->ld[i];
  }
}
#ifdef MPI
int dhist_process_l(dhist_t *d, int nfrm, MPI_File in)
{
  float *l;
  int st, i, err;

  l = (float *) calloc(d->fn, sizeof(float));
  if (!l)
    return -1;

  st = 3 * sizeof(int) + 4 * d->nfun * sizeof(double) + d->f0 * sizeof(float);
  for (i = 0; i < nfrm; i++) {
    err = MPI_File_read_at(in, st, l, d->fn, MPI_FLOAT, MPI_STATUS_IGNORE);
    if (err)
      goto err2;
    dhist_update(l, d);
    st += d->nfun * sizeof(float);
  }
  free(l);
  return 0;
 err2:
  set_mpi_error(err);  
 err1:
  free(l);
  return -1;
}
#else
int dhist_process_l(dhist_t *d, int nfrm, FILE *in)
{
  float *l;
  int st, i;

  l = (float *) calloc(d->fn, sizeof(float));
  if (!l)
    return -1;

  st = d->nfun - d->fn;
  if (nfrm) {
    if (fseek(in, d->f0 * sizeof(float), SEEK_CUR))
      goto err1;
    if (fread(l, sizeof(float), d->fn, in) != d->fn)
      goto err1;
    dhist_update(l, d);
  }
  for (i = 1; i < nfrm; i++) {
    if (fseek(in, st * sizeof(float), SEEK_CUR))
      goto err1;
    if (fread(l, sizeof(float), d->fn, in) != d->fn)
      goto err1;
    dhist_update(l, d);
  }
  free(l);
  return 0;
 err1:
  free(l);
  return -1;
}
#endif

#ifdef MPI
int dhist_write_hdr(dhist_t *d, MPI_File f)
{
  int err;
  err = MPI_File_write_at(f, 0, &d->dr, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
  if (err)
    goto mpi_err1;
  err = MPI_File_write_at(f, sizeof(double), &d->np, 4, MPI_INT, MPI_STATUS_IGNORE);
  if (err)
    goto mpi_err1;

  return 0;
 mpi_err1:
  set_mpi_error(err);
  return -1;
}
#else
int dhist_write_hdr(dhist_t *d, FILE *f)
{
  if (fwrite(&d->dr, sizeof(double), 1, f) != 1)
    return -1;
  if (fwrite(&d->np, sizeof(int), 4, f) != 4)
    return -1;
  return 0;
}
#endif

#ifdef MPI
int dhist_read_hdr(dhist_t *d, int n, MPI_File f)
{
  int m, l, i, err;
  err = MPI_File_read_all(f, &d->dr, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
  if (err)
    goto err1;
  err = MPI_File_read_all(f, &d->np, 4, MPI_INT, MPI_STATUS_IGNORE);
  if (err)
    goto err1;

  m = 2 * d->nfun;
  d->lm = (int *) calloc(m, sizeof(int));
  if (!d->lm)
    goto err2;
  d->ld = d->lm + d->nfun;


  err = MPI_File_read_all(f, d->lm, m, MPI_INT, MPI_STATUS_IGNORE);
  if (err)
    goto err1;
  
  m = 0;
  for (i = 0; i < n; i++)
    m += d->ld[i];
  l = m;
  for (; i < d->nfun; i++) {
    m += d->ld[i] - d->ld[i - n];
    if (m > l) l = m;
  }
   
  d->hst = (unsigned *) calloc(l, sizeof(unsigned));
  if (!d->hst) {
    free(d->lm);
    goto err2;
  }
  return 0;
 err1:
  set_mpi_error(err); 
  return -1;
 err2:
  set_std_error();
  return -1;
}
#else
int dhist_read_hdr(dhist_t *d, int n, FILE *f)
{
  int m, l, i;
  if (fread(&d->dr, sizeof(double), 1, f) != 1)
    goto err1;
  if (fread(&d->np, sizeof(int), 4, f) != 4)
    goto err1;

  m = 2 * d->nfun;
  d->lm = (int *) calloc(m, sizeof(int));
  if (!d->lm)
    goto err1;
  d->ld = d->lm + d->nfun;
  if (fread(d->lm, sizeof(int), m, f) != m) {
    free(d->lm);
    goto err1;
  }
   
  m = 0;
  for (i = 0; i < n; i++)
    m += d->ld[i];
  l = m;
  for (; i < d->nfun; i++) {
    m += d->ld[i] - d->ld[i - n];
    if (m > l) l = m;
  }
   
  d->hst = (unsigned *) calloc(l, sizeof(unsigned));
  if (!d->hst) {
    free(d->lm);
    goto err1;
  }
  return 0;
 err1:
  set_std_error();
  return -1;
}
#endif

#ifdef MPI
int dhist_write_hist(dhist_t *d, MPI_File f)
{
  int err, st;

  st = sizeof(double) + (4 + d->f0) * sizeof(int);
  err = MPI_File_write_at_all(f, st, d->lm, d->fn, MPI_INT, MPI_STATUS_IGNORE);
  if (err)
    goto mpi_err2;

  st += d->nfun * sizeof(int);
  err = MPI_File_write_at_all(f, st, d->ld, d->fn, MPI_INT, MPI_STATUS_IGNORE);
  if (err)
    goto mpi_err2;

  st += (d->nfun - d->f0) * sizeof(int) + d->hst0 * sizeof(unsigned);
  err = MPI_File_write_at_all(f, st, d->hst, d->hstn, MPI_UNSIGNED, MPI_STATUS_IGNORE);
  if (err)
    goto mpi_err2;

  return 0;
 mpi_err2:
  set_mpi_error(err); 
  return -1;
}
#else
int dhist_write_hist(dhist_t *d, FILE *f) 
{
  int m;
  m = 2 * d->fn;
  if (fwrite(d->lm, sizeof(int), m, f) != m)
    return -1;

  if (fwrite(d->hst, sizeof(unsigned), d->hstn, f) != d->hstn)
    return-1;

  return 0;
}
#endif
