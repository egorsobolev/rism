#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <math.h>
#include <float.h>
#include <cblas.h>

#include "grid.h"
#include "water.h"
#include "mol.h"
#include "avgw.h"
#include "eqoz.h"
#include "closure.h"
#include "nr.h"
#include "tdyn.h"

#include "tstutil.h"

#include "opt.h"

int awrism_mgrid(eq_t *eq, grid_param_t *gp, mol_t *m)
{
  int i, j, l, u, exitcode, nfun, np, maxit, incu;
  water_t w;
  grid_t g;
  int nn[3], n;
  rismaw_t *rism = (rismaw_t *) eq->p;

  double *tuv, *cuv, rtol, etol;

  nfun = m->natom * (m->natom - 1) / 2;


  for (i = 0; i < gp->n; i++) {
    if (water_read(gp->f[i].solv_file, &w)) {
      printf("Error reading file of water data '%s'\n", gp->f[i].solv_file);
      return 3;
    }
    printf("3.%d. SOLVENT\n", i+1);
    printf("   Natom     Nfun    Ngrid        T,K   rho,A^-3\n");
    printf("%8d %8d %8d %10g %10g\n", w.natom, w.nfun, w.ngrid, w.m.t / 1.9872e-3, w.m.rho);
 
    rism->natv = w.natom;
    rism->nfun = rism->natv * rism->natu;
    rism->v.xvv = w.xvv;
    rism->v.symc = w.symc;
    rism->v.n = w.n;
    rism->v.ngrid = w.ngrid;
    rism->v.t = w.m.t;
    rism->v.rho = w.m.rho;

    if (ginit(w.ngrid, w.dr, &g)) {
    }
    printf("4.%d. GRID\n", i+1);
    printf(" Ngrid      dr, A    dk, 1/A       L, A\n");
    printf("%6d %10g %10g %10g\n", g.n, g.dr, g.dk, g.dr * g.n);

    exitcode = avgw_read(gp->f[i].avgw_file, &rism->wuu, w.ngrid - 1, rism->reduc);
    if (exitcode) {
      printf("Error reading file of average omega matrix '%s'\n", gp->f[i].avgw_file);
      printf("code = %d\n", exitcode);
      exitcode = 4;
      goto err1;
    }
    printf("5.%d. AVERAGE OMEGA\n", i+1);
    printf(" Ngrid\n");
    printf("%6d\n", rism->wuu.np);
    if (((float) rism->wuu.dr != (float) g.dr) || (rism->wuu.np != g.n)) {
      printf("Grids of solvent functions and intramolecualar functions does not correspond\n");
      printf("ngrid=%d dr=%g\n", rism->wuu.np, rism->wuu.dr);
      exitcode = 5;
      goto err2;
    }
    if (rism->wuu.nfun != nfun) {
      printf("Nfun in the average matrix does not correspond to the Natom in molecule\n");

      printf("Nfun = %d\n", rism->wuu.nfun);
      exitcode = 6;
      goto err2;
    }
    
    if (poten_mk(&g, &w, m, rism->ngalpha, &rism->puv)) {
      printf("Error while initialize potential functions\n");
      exitcode = 5;
      goto err2;
    }

    if (grid_eq_mk(g.n - 1, g.dr, &rism->ge)) {
      printf("Error while initialize equation FFT plan\n");
      exitcode = 6;
      goto err3;
    }
    if (grid_jac_mk(g.n / rism->reduc - 1, g.dr, &rism->gj)) {
      printf("Error while initialize Jacobi matrix FFT plan\n");
      exitcode = 7;
      goto err4;
    }

    eq->nZ = rism->nfun * rism->ge.n;
    eq->nJx = rism->nfun * rism->gj.n;

    rism->dcdt = (float *) calloc(eq->nJx, sizeof(float));
    if (!rism->dcdt) {
      printf("Error while allocating dcdt\n");
      exitcode = 8;
      goto err5;
    }

    tuv = (double *) calloc(2 * eq->nZ, sizeof(double));
    if (!tuv) {
      printf("Error while allocationg tuv\n");
      exitcode = 9;
      goto err6;
    }
    cuv = tuv + eq->nZ;

    incu = rism->natv * rism->ge.n;

    l = 0;
    for (u = 0; u < rism->natu; u++) {
      j = rism->puv.atyp[u] * incu;
      for (i = 0; i < incu; i++) {
	tuv[l] = -rism->puv.asympr[j + i];
	l++;
      }
    }

    nn[0] = g.n - 1;
    nn[1] = w.natom;
    nn[2] = m->natom;
    write_dmtx("tuv.bin", 3, nn, tuv);


    maxit = eq->nrprm.maxi;
    rtol = eq->nrprm.rtol;
    etol = eq->nrprm.etol;

    exitcode = nr(eq, tuv, &rtol, &etol, &maxit);
    if (exitcode) {
      printf("NR error: %d\n", exitcode);
      goto err7;
    }

    rism->closure_c(rism, tuv, cuv);
    cblas_daxpy(eq->nZ, 1.0, cuv, 1, tuv, 1);

    printf("Etot = %.2f (kcal/mol)\n", etot(rism, tuv, cuv));
    printf("mu[SC/HNC] = %.2f (kcal/mol)\n", musc_hnc(rism, tuv, cuv));
    printf("mu[SC/PLHNC] = %.2f (kcal/mol)\n", musc_plhnc(rism, tuv, cuv));
    printf("mu[GF] = %.2f (kcal/mol)\n", mugf(rism, tuv, cuv));

    /*
    l = 0;
    for (u = 0; u < rism->natu; u++) {
      j = rism->puv.atyp[u] * incu;
      for (i = 0; i < incu; i++) {
	tuv[l] += rism->puv.asympr[j + i];
	l++;
      }
    }
    nn[0] = g.n - 1;
    nn[1] = w.natom;
    nn[2] = m->natom;
    write_dmtx("tuv.bin", 3, nn, tuv);
    */
    free(tuv);
    free(rism->dcdt);
    grid_jac_del(&rism->gj);
    grid_eq_del(&rism->ge);
    poten_del(&rism->puv);
    water_del(&w);
    avgw_del(&rism->wuu);
  }
  return 0;

 err7:
  free(tuv);
 err6:
  free(rism->dcdt);
 err5:
  grid_jac_del(&rism->gj);
 err4:
  grid_eq_del(&rism->ge);
 err3:
  poten_del(&rism->puv);
 err2:
  avgw_del(&rism->wuu);
 err1:
  water_del(&w);
  return exitcode;
}

int main(int narg, char **argv)
{
  rism_opt_t opt;
  void **argtable;
  int nerrors;
  int exitcode;

  grid_param_t gp;
  mol_t m;

  int i;
  char *c;
  char *solv_file, *avgw_file;

  eq_t eq;
  rismaw_t rism;

  argtable = argtable_mk(&opt);

  if (!argtable) {
    fprintf(stderr, "%s: insufficient memory\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  nerrors = arg_parse(narg, argv, argtable);
  if (opt.help->count > 0) {
    printf("Usage: %s", argv[0]);
    arg_print_syntax(stdout, argtable,"\n");
    arg_print_glossary(stdout, argtable,"  %-26s %s\n");
    exitcode = EXIT_SUCCESS;
    goto err1;
  }
  /* special case: '--version' takes precedence error reporting */
  if (opt.ver->count > 0){
    printf("awrism 1.0 Feb 20, 2014\n");
    printf("Author: Egor Sobolev\n");
    printf("Copyright (C) 2004-2013 Institute of Mathematical Problems of Biology RAS\n");
    exitcode = EXIT_SUCCESS;
    goto err1;
  }
  /* If the parser returned any errors then display them and exit */
  if (nerrors > 0) {
    /* Display the error details contained in the arg_end struct.*/
    arg_print_errors(stderr, opt.end, argv[0]);
    fprintf(stderr, "Try '%s --help' for more information.\n", argv[0]);
    exitcode = 1;
    goto err1;
  }
  printf("1. GRID PARAMETERS\n");
  printf(" # solvent file                        average matrix file\n");
  gp.n = opt.grid->count;
  for (i = 0; i < gp.n; i++) {
    c = strchr(opt.grid->sval[i], GRID_SEP);
    if (!c) {
      fprintf(stderr, "Invalid -g option.\n");
      fprintf(stderr, "Try '%s --help' for more information.\n", argv[0]);
      exitcode = 1;
      goto err1;
    }
    gp.f[i].solv_file = c + 1;
    gp.f[i].avgw_file = opt.grid->sval[i];
    *c = '\0';
    gp.o[i] = i;

    printf("%2d %-35.35s %-35.35s\n", i + 1, gp.f[i].avgw_file, gp.f[i].solv_file);
  }

  if (!strcmp(opt.closure->sval[0], "HNC")) {
    printf("Closure: HNC\n");
    rism.closure = &hnc;
    rism.closure_c = &hnc_c;
  } else if (!strcmp(opt.closure->sval[0], "PLHNC")) {
    printf("Closure: PLHNC\n");
    rism.closure = &plhnc;
    rism.closure_c = &plhnc_c;
  } else {
    fprintf(stderr, "Invalid -c option.\n");
    fprintf(stderr, "Try '%s --help' for more information.\n", argv[0]);
    goto err1;
  }
  if (opt.reduc->ival[0] < 1) {
    fprintf(stderr, "Invalid --reduc option. It must be equal to or greater then 1\n");
    fprintf(stderr, "Try '%s --help' for more information.\n", argv[0]);
    goto err1;
  }

  if (mol_read(opt.mol->filename[0], &m)) {
    printf("Error reading file of molecule '%s'\n", opt.mol->filename[0]);
    exitcode = 2;
    goto err1;
  }

  printf("2. MOLECULE\n");
  printf("   Natom    Nuniq  sum(q[i])\n");
  printf("%8d %8d %10f\n", m.natom, m.ntype, m.qtot);

  rism.natu = m.natom;
  rism.reduc = opt.reduc->ival[0];
  rism.ngalpha = opt.ngalpha->dval[0];

  eq.p = &rism;
  eq.Jx = &rismaw_Jx;
  eq.Z = &rismaw_eq;
  eq.getb = &rismaw_getx;
  eq.putx = &rismaw_putx;
  
  eq.nrprm.maxi = opt.nrmaxi->ival[0] - 1;
  eq.nrprm.etol = opt.etol->dval[0];
  eq.nrprm.rtol = opt.rtol->dval[0];
  eq.nrprm.lmaxi = opt.lmaxi->ival[0];
  eq.nrprm.narm = opt.narm->ival[0];

  awrism_mgrid(&eq, &gp, &m);

  
  exitcode = EXIT_SUCCESS;
 err3:
 err2:
  mol_del(&m);
 err1:
  argtable_del(argtable);
  exit(exitcode);
}
