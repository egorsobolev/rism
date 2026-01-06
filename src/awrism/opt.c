#include "opt.h"

#include <stdlib.h>

#define NOPTION 15

void **argtable_mk(rism_opt_t *opt)
{
  int i = 0;
  void **t;
  t = calloc(NOPTION, sizeof(void*));
  if (!t)
    return NULL;
  t[i++] = opt->mol = arg_file1("p", "mol", "<mol>", "molecule parameter file");
  t[i++] = opt->grid = arg_strn("g", "grid", "<avgw>:<solv>", 1, GRID_MAX, "grid specific parameters:");
  t[i++] = arg_rem(NULL, "<avgw> is a file that sets average matrix Omega");
  t[i++] = arg_rem(NULL, "<solv> is a file that sets the solvent data");
  t[i++] = opt->closure = arg_str0("c", NULL, "<closure>", "closure equation: HNC, PLHNC (default: PLHNC)");
  t[i++] = opt->reduc = arg_int0(NULL, "reduc", "<int>", "grid reduction times of linear solver (default: 4)");
  t[i++] = opt->ngalpha = arg_dbl0(NULL, "ngalpha", "<float>", "renormalization constant of Coulomb potential (default: 1.08)");
  t[i++] = opt->rtol = arg_dbl0(NULL, "rtol", "<float>", "residual RMSE tolerance (default: 1e-7)");
  t[i++] = opt->etol = arg_dbl0(NULL, "etol", "<float>", "internal energy error tolerance, kcal/mol (default: 0.2)");
  t[i++] = opt->nrmaxi = arg_int0(NULL, "nrmaxi", "<int>", "maximum number of NR iterations (default: 300)");
  t[i++] = opt->lmaxi = arg_int0(NULL, "lmaxi", "<int>", "maximum number of linear solver iterations (default: 5000)");
  t[i++] = opt->narm = arg_int0(NULL, "narm", "<int>", "maximum number of Armijo' step-size procedure iterations (default: 5)");
  t[i++] = opt->help = arg_lit0(NULL, "help", "print this help and exit");
  t[i++] = opt->prefix = arg_file1(NULL, NULL, "<prefix>", "prefix of output files");
  t[i++] = opt->end = arg_end(20);

  if (arg_nullcheck(t) != 0) {
    arg_freetable(t, NOPTION);
    free(t);
    return NULL;
  }

  opt->closure->sval[0] = "PLHNC";
  opt->reduc->ival[0] = 4;
  opt->ngalpha->dval[0] = 1.08;
  opt->rtol->dval[0] = 1e-7;
  opt->etol->dval[0] = 0.2;
  opt->nrmaxi->ival[0] = 300;
  opt->lmaxi->ival[0] = 5000;
  opt->narm->ival[0] = 5;

  return t;
}

void argtable_del(void **t)
{
  arg_freetable(t, NOPTION);
  free(t);
}
