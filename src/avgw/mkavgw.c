#include <stdlib.h>
#include <stdio.h>

#include <rism/dhist.h>
#include <rism/grid.h>
#include <rism/interp.h>

#include <argtable2.h>

#include <time.h>
#define walltime()      ((double) clock() / CLOCKS_PER_SEC)

#define GRID_SHAPES_MAX 16
#include "avgw.h"

int main(int narg, char **argv)
{
  struct arg_str *opt_shape;
  struct arg_lit *opt_help, *opt_ver;
  struct arg_file *opt_in, *opt_out;
  struct arg_end *end;
  void *argtable[] = {
    opt_shape = arg_strn("s", "shape", "<tr>:<np>:<dr>", 1, GRID_SHAPES_MAX, "final grid shape:"),
    arg_rem(NULL, "<tr> is cut off a fraction of the tail function"),
    arg_rem(NULL, "<np> is size in points"),
    arg_rem(NULL, "<dr> is spatial grid step, A"),
    opt_help = arg_lit0(NULL, "help", "print this help and exit"),
    opt_ver = arg_lit0(NULL, "version", "print version information and exit"),
    opt_in = arg_file1(NULL, NULL, "<h-file>", "input histograms file"),
    opt_out = arg_file1(NULL, NULL, "<aw-file>", "output aw-files prefix"),
    end = arg_end(20),
  };
  int nerrors;
  FILE *in;
  double tic, toc;
  int exitcode;
  avgw_shapes_t s;
  avgw_mtx_t W[GRID_SHAPES_MAX];
  avgw_cutparam_t awcut[GRID_SHAPES_MAX];
  sgrid_t g;
  dhist_t d;
  int i, j;
  avgw_func_t f;
  float *fcut;

  tic = walltime();
  if (arg_nullcheck(argtable) != 0) {
    /* NULL entries were detected, some allocations must have failed */
    fprintf(stderr, "%s: insufficient memory\n", argv[0]);
    exitcode = EXIT_FAILURE;
    goto err1;
  }

  nerrors = arg_parse(narg, argv, argtable);
  if (opt_help->count > 0) {
    printf("Usage: %s", argv[0]);
    arg_print_syntax(stdout, argtable,"\n");
    arg_print_glossary(stdout, argtable,"  %-26s %s\n");
    exitcode = EXIT_SUCCESS;
    goto err1;
  }
  /* special case: '--version' takes precedence error reporting */
  if (opt_ver->count > 0){
    printf("avgw 1.0 Aug 11, 2013\n");
    printf("Author: Egor Sobolev\n");
    printf("Copyright (C) 2004-2013 Institute of Mathematical Problems of Biology RAS\n");
    exitcode = EXIT_SUCCESS;
    goto err1;
  }
  exitcode = EXIT_FAILURE;
  /* If the parser returned any errors then display them and exit */
  if (nerrors > 0) {
    /* Display the error details contained in the arg_end struct.*/
    arg_print_errors(stderr, end, argv[0]);
    fprintf(stderr, "Try '%s --help' for more information.\n", argv[0]);
    goto err1;
  }


  in = fopen(opt_in->filename[0], "rb");
  if (!in) {
    perror("fopen");
    goto err1;
  }
  if (dhist_read_hdr(&d, in)) {
    perror("dhist_read_hdr");
    goto err2;
  }
  
  if (sgrid_init(d.np, d.dr, &g)) {
    perror("grid_init");
    goto err3;
  }

  s.interp = 0;
  s.n = opt_shape->count;
  for (i = 0; i < s.n; i++) {
    if (avgw_shape_parse(&g, opt_shape->sval[i], &s.p[i])) {
      fprintf(stderr, "%s: invalid value of option ", argv[0]);
      arg_print_option(stderr, "s", "shape", "<tr>:<np>:<dr>", "\n");
      fprintf(stderr, "Try '%s --help' for more information.\n", argv[0]);
      goto err4;
    }
    s.interp |= s.p[i].interp;
    j = i;
    while (j > 0 && (s.p[i].trun < s.p[s.o[j - 1]].trun || s.p[i].trun == s.p[s.o[j - 1]].trun && s.p[i].l > s.p[s.o[j - 1]].l)) {
      s.o[j] = s.o[j - 1];
      j--;
    }
    s.o[j] = i;
  }
  printf("Hist: np = %d, dR = %g (A), L = %g (A)%s\n", g.np, g.dr, g.np * g.dr, s.interp ? ", int": "");
  avgw_shapes_print(&s);

  if (avgw_func_init(&g, &f)) {
    perror("avgw_func_init");
    goto err4;
  }
  fcut = (float *) calloc(g.np, sizeof(float));
  if (!fcut) {
    perror("calloc");
    goto err5;
  }
  if (avgw_mtx_init(s.n, s.p, d.nfun, W)) {
    perror("avgw_mtx_init");
    goto err6;
  }
  if (avgw_shapes_open(opt_out->filename[0], s.n, s.p, s.f)) {
    perror("avgw_shapes_open");
    goto err7;
  }
  for (i = 0; i < s.n; i++) {
    if (fseek(s.f[i], (W[i].nfun + 2) * sizeof(int) + sizeof(float), SEEK_SET)) {
      perror("fseek");
      goto err8;
    }
  }

  for (j = 0; j < d.nfun; j++) {
    if (fread(d.hst, sizeof(unsigned), d.ld[j], in) != d.ld[j]) {
      perror("fread");
      goto err8;
    }
    avgw_hist2aw(&g, d.ld[j], d.lm[j], d.nfrm, d.hst, &f);
    avgw_itail(&g, &f, &s, awcut);
    if (s.interp)
      sspline_uni(f.np, f.s, 0.0, 0.0, f.s2, f.u);

    for (i = 0; i < s.n; i++) {
      W[i].Icut[j] = awcut[i].Icut;
      W[i].npcut[j] = awcut[i].npcut;
      W[i].nz[j] = awcut[i].nz;
      avgw_reshape(&g, &f, &s.p[i], &awcut[i], fcut);
      if (fwrite(fcut, sizeof(float), awcut[i].npcut, s.f[i]) != awcut[i].npcut) {
	perror("fwrite");
	goto err8;
      }
    }
  }
  for (i = 0; i < s.n; i++) {
    if (avgw_write_info(&W[i], s.f[i])) {
      perror("avgw_write_hdr");
      goto err8;
    }
    if (fseek(s.f[i], 0, SEEK_SET)) {
      perror("fseek");
      goto err8;
    }
    if (avgw_write_hdr(&W[i], s.f[i])) {
      perror("avgw_write_hdr");
      goto err8;
    }
  }

  toc = walltime();
  printf("Time: %lg s\n", toc - tic);
  exitcode = EXIT_SUCCESS;
 err8:
  avgw_shapes_close(s.n, s.f);
 err7:
  avgw_mtx_free(s.n, W);
 err6:
  free(fcut);
 err5:
  avgw_func_free(&f);
 err4:
  sgrid_free(&g);
 err3:
  dhist_free(&d);
 err2:
  fclose(in);
 err1:
  arg_freetable(argtable, sizeof(argtable)/sizeof(argtable[0]));
  exit(exitcode);
}
