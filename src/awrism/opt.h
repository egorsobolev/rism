#ifndef __RISM_AWRISM_OPT_H
#define __RISM_AWRISM_OPT_H

#include <argtable2.h>

struct RISM_OPTIONS
{
  struct arg_str *grid, *closure;
  struct arg_lit *help, *ver;
  struct arg_file *prefix, *mol;
  struct arg_dbl *ngalpha, *rtol, *etol;
  struct arg_int *reduc, *nrmaxi, *lmaxi;
  struct arg_int *narm;
  struct arg_end *end;
};
typedef struct RISM_OPTIONS rism_opt_t;

#define GRID_MAX 16
#define GRID_SEP ':'

struct GridFiles
{
  const char *avgw_file;
  const char *solv_file;
};

struct GridParam
{
  int n;
  int o[GRID_MAX];
  struct GridFiles f[GRID_MAX];
};
typedef struct GridParam grid_param_t;


void **argtable_mk(rism_opt_t *);
void argtable_del(void **);

#endif /* __RISM_AWRISM_OPT_H */
