#include <errno.h>
#include <string.h>
#include <stdio.h>

static int errcode;
static char *errstr;

#ifdef MPI
#include <mpi/mpi.h>
static char mpierr[MPI_MAX_ERROR_STRING];

void set_mpi_error(int e)
{
  int l;
  errcode = e;
  MPI_Error_string(e, mpierr, &l);
  errstr = mpierr;
}
#endif

void set_std_error()
{
  errcode = errno;
  errstr = strerror(errno);
}

void print_error(const char *m)
{
  fprintf(stderr, "%s: %s\n", m, errstr);
}
