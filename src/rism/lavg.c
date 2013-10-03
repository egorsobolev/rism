#include <stdlib.h>
#include <memory.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include <rism/lavg.h>

#ifdef MPI
#include <mpi.h>
#endif

int lavg_init(int natm, lavg_t *s)
{
  int i;
  s->natm = natm;
  s->nfun = natm * (natm - 1) / 2;
  s->nfrm = 0;
  s->mn = (double *) calloc(4 * s->nfun, sizeof(double));
  if (!s->mn) {
    set_std_error();
    return -1;
  }
  s->mx = s->mn + s->nfun;
  s->mean = s->mx + s->nfun;
  s->std = s->mean + s->nfun;
  for (i = 0; i < s->nfun; i++)
    s->mn[i] = DBL_MAX;
  return 0;
}

void lavg_free(lavg_t *s)
{
  free(s->mn);
}

#ifdef MPI
int lavg_writehdr(const lavg_t *s, MPI_File f)
{
  int n, err;
  n = 4 * s->nfun;


  err = MPI_File_write_at(f, 0 , s, 3, MPI_INT, MPI_STATUS_IGNORE);
  if (err) goto mpi_err1;
  err = MPI_File_write_at(f, 3 * sizeof(int), s->mn, n, MPI_DOUBLE, MPI_STATUS_IGNORE);
  if (err) goto mpi_err1;

  return 0;
 mpi_err1:
  set_mpi_error(err);
  return -1;
}
#else
int lavg_writehdr(const lavg_t *s, FILE *f)
{
  int n;
  n = 4 * s->nfun;
  fseek(f, 0L, SEEK_SET);
  if (fwrite(s, sizeof(int), 3, f) != 3)
    goto err1;
  if (fwrite(s->mn, sizeof(double), n, f) != n)
    goto err1;

  return 0;
 err1:
  set_std_error();
  return -1;
}
#endif

#ifdef MPI
int lavg_readhdr(lavg_t *s, MPI_File f)
{
  int n, err;
  err = MPI_File_read_all(f, s, 3, MPI_INT, MPI_STATUS_IGNORE);
  if (err)
    goto mpi_err1;

  n = 4 * s->nfun;
  s->mn = (double *) calloc(n, sizeof(double));
  if (!s->mn)
    return -1;
  s->mx = s->mn + s->nfun;
  s->mean = s->mx + s->nfun;
  s->std = s->mean + s->nfun;

  err = MPI_File_read_all(f, s->mn, n, MPI_DOUBLE, MPI_STATUS_IGNORE);
  if (err)
    goto mpi_err1;

  return 0;
 mpi_err1:
  set_mpi_error(err);  
  return -1;
}
#else
int lavg_readhdr(lavg_t *s, FILE *f)
{
  int n;
  if (fread(s, sizeof(int), 3, f) != 3)
    return -1;
  n = 4 * s->nfun;
  s->mn = (double *) calloc(n, sizeof(double));
  if (!s->mn)
    return -1;
  s->mx = s->mn + s->nfun;
  s->mean = s->mx + s->nfun;
  s->std = s->mean + s->nfun;
  if (fread(s->mn, sizeof(double), n, f) != n)
    return -1;
  return 0;
}
#endif

void lavg_update(const float *x, float *l, lavg_t *s)
{
  int k, i, j, n;
  double dx, dy, dz, a;
  n = 3 * s->natm;
  k = 0;
  for (i = 3; i < n; i += 3) {
    for (j = 0; j < i; j += 3) {
      dx = x[i] - x[j];
      dy = x[i + 1] - x[j + 1];
      dz = x[i + 2] - x[j + 2];
      a = dx * dx + dy * dy + dz * dz;
      s->std[k] += a;
      a = sqrt(a);
      l[k] = a;
      s->mean[k] += a;
      if (a < s->mn[k])
	s->mn[k] = a;
      if (a > s->mx[k])
	s->mx[k] = a;
      k++;
    }
  }
  s->nfrm++;
}

int lavg_finish(lavg_t *s)
{
  int i, nfrm;
  double m;

#ifdef MPI
  int rank, err;
  MPI_Allreduce(MPI_IN_PLACE, &s->nfrm, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (!s->nfrm)
    return 0;

#ifdef MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank) {
    err = MPI_Reduce(s->mn, NULL, s->nfun, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    if (err) goto mpi_err1;
    err = MPI_Reduce(s->mx, NULL, s->nfun, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (err) goto mpi_err1;
    err = MPI_Reduce(s->mean, NULL, 2 * s->nfun, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (err) goto mpi_err1;
    return 0;
  } else {
    err = MPI_Reduce(MPI_IN_PLACE, s->mn, s->nfun, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    if (err) goto mpi_err1;
    err = MPI_Reduce(MPI_IN_PLACE, s->mx, s->nfun, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (err) goto mpi_err1;
    err = MPI_Reduce(MPI_IN_PLACE, s->mean, 2 * s->nfun, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (err) goto mpi_err1;
  }
#endif

  if (s->nfrm == 1) {
    memset(s->std, 0, s->nfun * sizeof(double));
    return 0;
  }

  m = 1.0 / (s->nfrm - 1);
  for (i = 0; i < s->nfun; i++) {
    s->mean[i] /= s->nfrm;
    s->std[i] = sqrt(m * (s->std[i] - s->nfrm * s->mean[i] * s->mean[i]));
  }
  return 0;
#ifdef MPI
 mpi_err1:
  set_mpi_error(err);
  return -1;
#endif
}
