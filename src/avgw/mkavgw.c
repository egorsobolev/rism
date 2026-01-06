#include <rism/dhist.h>
#include <rism/grid.h>
#include <rism/interp.h>
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

#define GRID_SHAPES_MAX 16
#include "avgw.h"

int ncycles(int nf, int ne, int np)
{
  int n, r, nsess;
  if (np > nf)
    return -1;
  n = ne * np;
  nsess = nf / n;
  r = nf - nsess * n;
  /* если остаток больше числа процессоров, то просто берем на один блок больше */
  /* иначе смотрим сколько блоков нужно, если размер блока уменьшить на 1 */
  if (r >= np) {
    nsess++;
  } else if (r) {
    nsess = nf / (n - np);
    r = nf - nsess * n;
    /* если числа блоков не хватает даже при использовании их исходного размера,
     * используем на один блок больше */
    if (r > 0)
      nsess++;
  }
  return nsess;
}

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
  double tic, toc;
  int exitcode;
  avgw_shapes_t s;
  avgw_mtx_t W[GRID_SHAPES_MAX];
  avgw_cutparam_t awcut[GRID_SHAPES_MAX];
  sgrid_t g;
  dhist_t d;
  int i, j, k;
  avgw_func_t f;
  /*  float *fcut; */
  avgw_outbuf_t cut[GRID_SHAPES_MAX];
  unsigned long pos;
  int nc, ns, cpos, bpos, blen, nb, r, hlen;
  int nproc, rank;

#ifdef MPI
  int err;
  MPI_File in;

  MPI_Init(&narg, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  if (!rank)
    printf("Run on %d processors\n", nproc);
#else
  FILE *in;

  nproc = 1;
  rank = 0;
#endif

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

  if (dhist_read_hdr(&d, AVGW_BUFSIZE, in)) {
    print_error("dhist_read_hdr");
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
  if (!rank) {
    printf("Hist: np = %d, dR = %g (A), L = %g (A)%s\n", g.np, g.dr, g.np * g.dr, s.interp ? ", int": "");
    avgw_shapes_print(&s);
  }

  if (avgw_func_init(&g, &f)) {
    perror("avgw_func_init");
    goto err4;
  }
  if (avgw_outbuf_init(s.n, s.p, cut)) {
    perror("avgw_outbuf_init");
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
#ifdef MPI
    err = MPI_File_seek_shared(s.f[i], (2 + d.nfun) * sizeof(int) + sizeof(float), MPI_SEEK_SET);
    if (err) {
      set_mpi_error(err);
      print_error("MPI_File_seek_shared:");
      goto err8;
    }
#else
    if (fseek(s.f[i], (2 + d.nfun) * sizeof(int) + sizeof(float), SEEK_SET)) {
      perror("fseek");
      goto err8;
    }
#endif
  }
#ifdef MPI
  err = MPI_File_seek_shared(in, sizeof(double) + (4 + 2 * d.nfun) * sizeof(int), MPI_SEEK_SET);
  if (err) {
    set_mpi_error(err);
    print_error("MPI_File_seek_shared:");
    goto err8;
  }
#endif
  nc = ncycles(d.nfun, AVGW_BUFSIZE, nproc);
  cpos = 0;
  for (j = 0; j < nc; j++) {
    ns = d.nfun / nc;
    r = d.nfun - ns * nc;
    ns += (j < r);

    blen = ns / nproc;
    r = ns - blen * nproc;
    bpos = cpos + blen * rank + (rank < r ? rank : r);
    blen += (rank < r);

    hlen = d.ld[bpos];
    for (i = 1; i < blen; i++)
      hlen += d.ld[bpos + i];

#ifdef MPI
    err = MPI_File_read_ordered(in, d.hst, hlen, MPI_UNSIGNED, MPI_STATUS_IGNORE);
    if (err) {
      set_mpi_error(err);
      print_error("MPI_File_read_ordered:");
      goto err8;
    }
#else
    if (fread(d.hst, sizeof(unsigned), hlen, in) != hlen) {
      perror("fread");
      goto err8;
    }
#endif
    r = 0;
    for (k = 0; k < blen; k++) {
      avgw_hist2aw(&g, d.ld[bpos], d.lm[bpos], d.nfrm, d.hst + r, &f);
      avgw_itail(&g, &f, &s, awcut);
      if (s.interp)
        sspline_uni(f.np, f.s, 0.0, 0.0, f.s2, f.u);

      for (i = 0; i < s.n; i++) {
        W[i].Icut[bpos] = awcut[i].Icut;
        W[i].npcut[bpos] = awcut[i].npcut;
        W[i].nz[bpos] = awcut[i].nz;

        avgw_reshape(&g, &f, &s.p[i], &awcut[i], cut[i].cur);
        cut[i].cur = cut[i].cur + awcut[i].npcut;
      }
      r += d.ld[bpos];
      bpos++;
    }

    for (i = 0; i < s.n; i++) {
      nb = cut[i].cur - cut[i].mem;
#ifdef MPI
      err = MPI_File_write_ordered(s.f[i], cut[i].mem, nb, MPI_FLOAT, MPI_STATUS_IGNORE);
      if (err) {
        set_mpi_error(err);
        print_error("MPI_File_write_ordered:");
        goto err8;
      }
#else
      if (fwrite(cut[i].mem, sizeof(float), nb, s.f[i]) != nb) {
        perror("fwrite");
        goto err8;
      }
#endif
      cut[i].cur = cut[i].mem;
    }
    cpos += ns;
  }

  for (i = 0; i < s.n; i++) {
#ifdef MPI
    err = MPI_File_seek_shared(s.f[i], 2 * sizeof(int) + sizeof(float), MPI_SEEK_SET);
    if (err) {
      set_mpi_error(err);
      print_error("MPI_File_seek_shared:");
      goto err8;
    }
#else
    if (fseek(s.f[i], 2 * sizeof(int) + sizeof(float), SEEK_SET)) {
      perror("fseek");
      goto err8;
    }
#endif
  }

#ifdef MPI
  cpos = 0;
  for (j = 0; j < nc; j++) {
    ns = d.nfun / nc;
    r = d.nfun - ns * nc;
    ns += (j < r);

    blen = ns / nproc;
    r = ns - blen * nproc;
    bpos = cpos + blen * rank + (rank < r ? rank : r);
    blen += (rank < r);

    for (i = 0; i < s.n; i++) {
      err = MPI_File_write_ordered(s.f[i], W[i].npcut + bpos, blen, MPI_INT, MPI_STATUS_IGNORE);
      if (err) {
        set_mpi_error(err);
        print_error("MPI_File_write_ordered:");
        goto err8;
      }
    }
    cpos += ns;
  }
#else
  for (i = 0; i < s.n; i++) {
    if (fwrite(W[i].npcut, sizeof(int), d.nfun, s.f[i]) != d.nfun) {
      perror("fwrite");
      goto err8;
    }
  }
#endif

  for (i = rank; i < s.n; i += nproc) {
    if (avgw_write_hdr(&W[i], s.f[i])) {
      print_error("avgw_write_hdr:");
      goto err8;
    }
  }

  /*
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
  */
  toc = walltime();
  if (!rank)
    printf("Time: %lg s\n", toc - tic);
  exitcode = EXIT_SUCCESS;
 err8:
  avgw_shapes_close(s.n, s.f);
 err7:
  avgw_mtx_free(W);
 err6:
  avgw_outbuf_free(cut);
  /* free(fcut); */
 err5:
  avgw_func_free(&f);
 err4:
  //sgrid_free(&g);
 err3:
  dhist_free(&d);
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
