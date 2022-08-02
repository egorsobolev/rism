#include "shapes.h"

#include <rism/grid.h>
#include <rism/strtools.h>
#include <rism/error.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

int avgw_shape_parse(const sgrid_t *g, const char *shape, gridparam_t *p)
{
  const char *b;
  char *e;
  double dt;
  b = shape;
  p->np = g->np;
  p->dr = g->dr;
  p->trun = strtod(b, &e);
  if (e == b)
    return -1;
  if (*e == ':') {
    b = e + 1;
    p->np = strtoul(b, &e, 10);
    if (e == b)
      p->np = g->np;
  }
  if (*e == ':') {
    b = e + 1;
    p->dr = strtod(b, &e);
    if (e == b)
      p->dr = g->dr;
  }
  if (*e)
    return -1;

  p->l = p->dr * p->np;
  p->dt = M_PI / p->l;
  p->interp = (fmodf(p->dt, g->dt) > FLT_EPSILON);

  return 0;
}

void avgw_shapes_print(avgw_shapes_t *s)
{
  int i;
  float trun;

  trun = NAN;
  printf("    trun     Np      dR(A)         L(A) int\n");
  for (i = 0; i < s->n; i++) {
    if (trun != s->p[s->o[i]].trun) {
      trun = s->p[s->o[i]].trun;
      printf("%8lg ", trun);
    } else
      printf("%9c", ' ');
    printf("%6d %10lg %12lg %s\n", s->p[s->o[i]].np, s->p[s->o[i]].dr, s->p[s->o[i]].l, s->p[s->o[i]].interp ? "yes" : "no");
  }
}
/*
int avgw_shapes_alloc(int n, const gridparam_t *p, void **b)
{
   int i, l;
   
   memset(b, 0, GRID_SHAPES_MAX * sizeof(float *));
   l = p[0].np;
   for (i = 1; i < n; i++)
     l += p[i].np;
   
   b[0] = malloc((l * sizeof(float) + sizeof(int)) * GRID_SHAPES_BUFSIZE);
   if (b[0])
     return -1;
   for (i = 1; i < n; i++)
      b[i] = (float *) ((int *) b[0] + GRID_SHAPE_BUFSIZE) + GRID_SHAPE_BUFSIZE * p[i-1].np;
   return 0;
}

void avgw_shapes_free(int n, void **b)
{
   free(b[0]);
}
*/
#ifdef MPI
int avgw_shapes_open(const char *prefix, int n, const gridparam_t *p, MPI_File *f)
{
  int i, l, j, err;
  char fn[4096];
  static char *ext = ".bin";

  l = strlen(prefix);
  if (l > 4060)
    l = 4060;
  for (i = 0; i < n; i++) {
    strncpy(fn, prefix, l);
    j = l;
    fn[j++] = '-';
    j += fractostr(fn + j, p[i].trun, FLT_EPSILON);
    fn[j++] = '-';
    j += sprintf(fn + j, "%d", p[i].np);
    fn[j++] = '-';
    j += fractostr(fn + j, p[i].dr, FLT_EPSILON);
    strcpy(fn + j, ext);

    err = MPI_File_open(MPI_COMM_WORLD, fn, MPI_MODE_WRONLY | MPI_MODE_CREATE , MPI_INFO_NULL, &f[i]);
    if (err)
      break;
    err = MPI_File_set_size(f[i], 0);
    if (err)
      break;
  }
  if (i == n)
    return 0;
  set_mpi_error(err);
  for (; i >=0; i--)
    MPI_File_close(&f[i]);
  return -1;
}

int avgw_shapes_close(int n, MPI_File *f)
{
  int i;
  for (i = 0; i < n; i++) {
    MPI_File_close(&f[i]);
  }
}
#else
int avgw_shapes_open(const char *prefix, int n, const gridparam_t *p, FILE **f)
{
  int i, l, j;
  char fn[4096];
  static char *ext = ".bin";

  l = strlen(prefix);
  if (l > 4060)
    l = 4060;
  for (i = 0; i < n; i++) {
    strncpy(fn, prefix, l);
    j = l;
    fn[j++] = '-';
    j += fractostr(fn + j, p[i].trun, FLT_EPSILON);
    fn[j++] = '-';
    j += sprintf(fn + j, "%d", p[i].np);
    fn[j++] = '-';
    j += fractostr(fn + j, p[i].dr, FLT_EPSILON);
    strcpy(fn + j, ext);

    f[i] = fopen(fn, "wb");
    if (!f[i])
      break;
  }
  if (i == n)
    return 0;
  for (; i >=0; i--)
    fclose(f[i]);
  return -1;
}

int avgw_shapes_close(int n, FILE **f)
{
  int i;
  for (i = 0; i < n; i++) {
    fclose(f[i]);
  }
}
#endif
