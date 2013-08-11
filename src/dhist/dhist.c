#include <stdlib.h>
#include <stdio.h>

#include <rism/lavg.h>
#include <rism/dhist.h>

#include <argtable2.h>

#include <time.h>
#define walltime()      ((double) clock() / CLOCKS_PER_SEC)

int main(int narg, char **argv)
{
  struct arg_int *opt_np;
  struct arg_dbl *opt_dr;
  struct arg_lit *opt_help, *opt_ver;
  struct arg_file *opt_in, *opt_out;
  struct arg_end *end;
  void *argtable[] = {
    opt_np = arg_int0("n", "points", "<np>", "grid size in points (default: 8192)"),
    opt_dr = arg_dbl0("s", "dr", "<dr>", "spatial grid step, A (default: 0.014)"),
    opt_help = arg_lit0(NULL, "help", "print this help and exit"),
    opt_ver = arg_lit0(NULL, "version", "print version information and exit"),
    opt_in = arg_file1(NULL, NULL, "<s-file>", "input atom-atom distance file"),
    opt_out = arg_file1(NULL, NULL, "<h-file>", "output histograms file"),
    end = arg_end(20),
  };
  int nerrors;
  FILE *in, *out;
  double tic, toc;
  int exitcode, i;
  lavg_t s;
  dhist_t d;
  float *l;

  tic = walltime();
  if (arg_nullcheck(argtable) != 0) {
    /* NULL entries were detected, some allocations must have failed */
    printf("%s: insufficient memory\n", argv[0]);
    exitcode = EXIT_FAILURE;
    goto err1;
  }

  nerrors = arg_parse(narg, argv, argtable);
  if (opt_help->count > 0) {
    printf("Usage: %s", argv[0]);
    arg_print_syntax(stdout, argtable,"\n");
    arg_print_glossary(stdout, argtable,"  %-20s %s\n");
    exitcode = EXIT_SUCCESS;
    goto err1;
  }
  opt_np->ival[0] = 8192;
  opt_dr->dval[0] = 0.014;
  /* special case: '--version' takes precedence error reporting */
  if (opt_ver->count > 0){
    printf("dhist 1.0 Aug 11, 2013\n");
    printf("Author: Egor Sobolev\n");
    printf("Copyright (C) 2004-2013 Institute of Mathematical Problems of Biology RAS\n");
    exitcode = EXIT_SUCCESS;
    goto err1;
  }
  exitcode = EXIT_FAILURE;
  /* If the parser returned any errors then display them and exit */
  if (nerrors > 0) {
    /* Display the error details contained in the arg_end struct.*/
    arg_print_errors(stdout, end, argv[0]);
    printf("Try '%s --help' for more information.\n", argv[0]);
    goto err1;
  }

  in = fopen(opt_in->filename[0], "rb");
  if (!in) {
    perror("fopen");
    goto err1;
  }
  if (lavg_readhdr(&s, in)) {
    perror("lavg_readhdr");
    goto err2;
  }

  if (dhist_init(opt_np->ival[0], opt_dr->dval[0], &s, &d)) {
    perror("dhist_init");
    goto err3;
  }

  l = (float *) malloc(s.nfun * sizeof(float));
  if (!l) {
    perror("malloc");
    goto err4;
  }

  for (i = 0; i < s.nfrm; i++) {
    if (fread(l, sizeof(float), s.nfun, in) != s.nfun) {
      perror("fread");
      goto err5;
    }
    dhist_update(l, &d);
  }

  out = fopen(opt_out->filename[0], "wb");
  if (!out) {
    perror("fopen");
    goto err5;
  }
  if (dhist_write(&d, out)) {
    perror("dhist_write");
    goto err6;
  }
  toc = walltime();
  printf("%10lf s, %d bytes, %d numbers\n", toc - tic, s.nfun * sizeof(float) + 2 * d.nfun * sizeof(int) + d.n * sizeof(unsigned), d.n);
  exitcode = EXIT_SUCCESS;

 err6:
  fclose(out);
 err5:
  free(l);
 err4:
  dhist_free(&d);
 err3:
  lavg_free(&s);
 err2:
  fclose(in);
 err1:
  arg_freetable(argtable, sizeof(argtable)/sizeof(argtable[0]));
  exit(exitcode);
}
