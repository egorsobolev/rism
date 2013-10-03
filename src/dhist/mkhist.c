#include <stdlib.h>
#include <stdio.h>

#include <rism/lavg.h>
#include <rism/dhist.h>

#include <argtable2.h>

#ifdef MPI
# include "mpi.h"
# define walltime()      MPI_Wtime()
#else
# include <time.h>
# define walltime()      ((double) clock() / CLOCKS_PER_SEC)
#endif

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
  double tic, toc;
  int exitcode;
  lavg_t s;
  dhist_t d;
  float *l;
  int rank, nproc;

#ifdef MPI
  int err;
  MPI_File in, out;

  MPI_Init(&narg, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  if (!rank)
    printf("Run on %d processors\n", nproc);
#else
  FILE *in, *out;

  nproc = 1;
  rank = 0;
#endif

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

#ifdef MPI
  err = MPI_File_open(MPI_COMM_WORLD, opt_in->filename[0], MPI_MODE_RDONLY , MPI_INFO_NULL, &in);
  if (err) {
    set_mpi_error(err);
    print_error("MPI_File_open:");
    goto err1;
  }  
#else
  in = fopen(opt_in->filename[0], "rb");
  if (!in) {
    perror("fopen");
    goto err1;
  }
#endif
  if (lavg_readhdr(&s, in)) {
    print_error("lavg_readhdr");
    goto err2;
  }

  if (dhist_init(opt_np->ival[0], opt_dr->dval[0], &s, nproc, rank, &d)) {
    perror("dhist_init");
    goto err3;
  }
  if (dhist_process_l(&d, s.nfrm, in)) {
    print_error("dhist_process_l");
    goto err3;
  }
#ifdef MPI
  err = MPI_File_open(MPI_COMM_WORLD, opt_out->filename[0], MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &out);
  if (!err)
    err = MPI_File_set_size(out, 0);
  if (err) {
    set_mpi_error(err);
    print_error("MPI_File_open:");
    goto err5;
  }
#else
  out = fopen(opt_out->filename[0], "wb");
  if (!out) {
    perror("fopen");
    goto err5;
  }
#endif
  if (!rank)
    if (dhist_write_hdr(&d, out)) {
      print_error("dhist_write_hdr");
      goto err6;
    }
  if (dhist_write_hist(&d, out)) {
    print_error("dhist_write_hist");
    goto err6;
  }
  toc = walltime();
  if (!rank)
    printf("%10lf s, %d bytes, %d numbers\n", toc - tic, s.nfun * sizeof(float) + 2 * d.nfun * sizeof(int) + d.n * sizeof(unsigned), d.n);
  exitcode = EXIT_SUCCESS;

 err6:
#ifdef MPI
  MPI_File_close(&out);
#else
  fclose(out);
#endif
 err5:
  dhist_free(&d);
 err3:
  lavg_free(&s);
 err2:
#ifdef MPI
  MPI_File_close(&in);
#else
  fclose(in);
#endif
 err1:
#ifdef MPI
  MPI_Finalize();
#endif
  arg_freetable(argtable, sizeof(argtable)/sizeof(argtable[0]));
  exit(exitcode);
}
