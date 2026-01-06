#include <rism/ambtraj.h>
#include <rism/lavg.h>
#include <rism/error.h>

#include <stdlib.h>
#include <stdio.h>

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
  struct arg_int *opt_natm, *opt_beg, *opt_step, *opt_fin;
  struct arg_lit *opt_ewald, *opt_help, *opt_ver;
  struct arg_file *opt_traj, *opt_out;
  struct arg_end *end;
  void *argtable[] = {
    opt_natm = arg_int1("a", "atoms", "<natm>", "number of atoms"),
    opt_ewald = arg_lit0("d", "ewald", "ewald box"),
    opt_beg = arg_int0("b", "start", "<f0>", "start from f0 frame (default: 1)"),
    opt_fin = arg_int0("e", "end", "<fn>", "finish at fn frame (default: until the end of file)"),
    opt_step = arg_int0("s", "step", "<m>", "process each m frame (default: 1)"),
    opt_help = arg_lit0(NULL, "help", "print this help and exit"),
    opt_traj = arg_file1(NULL, NULL, "<trajectory>", "AMBER trajectory file"),
    opt_out = arg_file1(NULL, NULL, "<output>", "Output file"),
    end = arg_end(20),
  };
  int nerrors;
  lavg_t s;
  float *x, *l;
  ambtraj_t traj;
  int f0, fn, nskip;
  double tic, toc;
  int exitcode;
  int rank, nproc;

#ifdef MPI
  MPI_File out;
  int err;
  long pos, blen;

  MPI_Init(&narg, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  if (!rank)
    printf("Run on %d processors\n", nproc);
#else
  FILE *out;

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
  opt_beg->ival[0] = 1;
  opt_fin->ival[0] = 0;
  opt_step->ival[0] = 1;

  nerrors = arg_parse(narg, argv, argtable);
  if (opt_help->count > 0) {
    printf("Usage: %s", argv[0]);
    arg_print_syntax(stdout, argtable,"\n");
    arg_print_glossary(stdout, argtable,"  %-20s %s\n");
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
  if (opt_beg->ival[0] < 1) {
    printf("invalid value of option -b, --start=<f0>\n");
    goto err1;
  }
  if (opt_fin->count > 0 && opt_fin->ival[0] < 1) {
    printf("invalid value of option -e, --end=<fn>\n");
    goto err1;
  }
  if (opt_step->ival[0] < 1) {
    printf("invalid value of option -s, --step=<m>\n");
    goto err1;
  }
  f0 = opt_beg->ival[0] + opt_step->ival[0] * rank - 1;
  fn = opt_fin->ival[0];
  nskip = opt_step->ival[0] * nproc - 1;
  if (fn > 0 && f0 >= fn) {
    printf("Nothing to do. Exit.\n");
    goto err1;
  }

  if (amb_open(opt_traj->filename[0], opt_natm->ival[0], opt_ewald->count > 0, &traj)) {
    perror("amb_open");
    goto err1;
  }
  if (!rank) {
    printf("Process trajectory\n");
    printf("%s\n", traj.title);
  }
  if (amb_ignore(&traj, f0) == -1) {
    perror("amb_ignore");
    goto err2;
  }
  if (!feof(traj.f)) {
    x = (float *) calloc(3 * opt_natm->ival[0], sizeof(float));
    if (!x) {
      perror("calloc");
      goto err2;
    }
    if (amb_read(&traj, x) == -1) {
      perror("amb_read");
      goto err6;
    }
  } else {
    printf("Trajectory has only %d frame(s).\n", traj.frm);
    printf("Nothing was done.\n");
    goto err3;
  }

  if (lavg_init(opt_natm->ival[0], &s)) {
    print_error("linit");
    goto err3;
  }
  l = (float *) calloc(s.nfun, sizeof(float));
  if (!l) {
    perror("calloc");
    goto err4;
  }

#ifdef MPI
  pos = 4 * s.nfun * sizeof(double) + 3 * sizeof(int) + s.nfun * sizeof(float) * rank;
  blen = nproc * s.nfun * sizeof(float);
  err = MPI_File_open(MPI_COMM_WORLD, opt_out->filename[0], MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &out);
  if (!err)
    err = MPI_File_set_size(out, 0);
  if (err) {
    set_mpi_error(err);
    print_error("MPI_File_open:");
    goto err5;
  }
  /*  MPI_File_set_view(out, 4 * s.nfun * sizeof(double) + 3 * sizeof(int) + s.nfun * sizeof(float) * rank, MPI_FLOAT, MPI_FLOAT, "native", MPI_INFO_NULL);*/
#else
  out = fopen(opt_out->filename[0], "wb");
  if (!out) {
    perror("fopen");
    goto err5;
  }
  fseek(out, 4 * s.nfun * sizeof(double) + 3 * sizeof(int), SEEK_SET);
#endif
  while (!feof(traj.f) && (fn <= 0 || traj.frm < fn)) {
    lavg_update(x, l, &s);
#ifdef MPI
    err = MPI_File_write_at(out, pos, l, s.nfun, MPI_FLOAT, MPI_STATUS_IGNORE);
    if (err) {
      set_mpi_error(err);
      print_error("lavg_writefrm");
      goto err6;
    }
    pos += blen;
#else
    if (fwrite(l, sizeof(float), s.nfun, out) != s.nfun) {
      perror("lavg_writefrm");
      goto err6;
    }
#endif
    if (!feof(traj.f) && (fn <= 0 || (traj.frm + nskip + 1) < fn)) {
      if (amb_ignore(&traj, nskip) == -1) {
        perror("amb_ignore");
        goto err6;
      }
      if (!feof(traj.f)) {
        if (amb_read(&traj, x) == -1) {
          perror("amb_read");
          goto err6;
        }
      }
    }
  }
  if (s.nfrm) {
    if (lavg_finish(&s)) {
      print_error("lavg_finish");
      goto err6;
    }
    if (!rank) {
      if (lavg_writehdr(&s, out)) {
        print_error("lavg_writehdr");
        goto err6;
      }
    }
    toc = walltime();
    if (!rank)
      printf("%10lf s, %d bytes, %d frames, %d pairs\n", toc - tic, (3 * s.natm + s.nfun) * sizeof(float) + 4 * s.nfun * sizeof(double), s.nfrm, s.nfun);
  } else {
    printf("Trajectory has only %d frame(s).\n", traj.frm);
    printf("Nothing was done.\n");
  }
  exitcode = EXIT_SUCCESS;
 err6:
#ifdef MPI
  MPI_File_close(&out);
#else
  fclose(out);
#endif
 err5:
  free(l);
 err4:
  lavg_free(&s);
 err3:
  free(x);
 err2:
  amb_close(&traj);
 err1:
#ifdef MPI
  MPI_Finalize();
#endif
  arg_freetable(argtable, sizeof(argtable)/sizeof(argtable[0]));
  exit(exitcode);
}
