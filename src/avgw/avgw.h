#ifndef __RISM_AVGW_AVGW_H
#define __RISM_AVGW_AVGW_H

#include "shapes.h"

#include <rism/grid.h>

#include <stdio.h>

#ifndef AVGW_MINZEROS
#define AVGW_MINZEROS 6
#endif

#ifndef AVGW_BUFSIZE
#define AVGW_BUFSIZE  4
#endif

struct AVGW_FUNC
{
	int np;
	int nz;
	double I;
	float *s;
	float *s2;
	float *u;
};
typedef struct AVGW_FUNC avgw_func_t;

struct AVGW_CUTPARAM
{
	int nz;
	int npcut;
	float Icut;
};
typedef struct AVGW_CUTPARAM avgw_cutparam_t;

struct AVGW_MTX
{
	int np;
	float dr;
	int nfun;
	int *npcut;
	int *nz;
	float *Icut;
};
typedef struct AVGW_MTX avgw_mtx_t;

struct AVGW_OUTBUF
{
	float *mem;
	float *cur;
};
typedef struct AVGW_OUTBUF avgw_outbuf_t;

int avgw_func_init(const grid_t *g, avgw_func_t *f);
void avgw_func_free(avgw_func_t *f);
int avgw_mtx_init(int n, gridparam_t *p, int nfun, avgw_mtx_t *W);
void avgw_mtx_free(avgw_mtx_t *W);
int avgw_outbuf_init(int n, const gridparam_t *p, avgw_outbuf_t *b);
int avgw_outbuf_free(avgw_outbuf_t *b);
#ifdef MPI
int avgw_write_hdr(avgw_mtx_t *W, MPI_File f);
#else
int avgw_write_hdr(avgw_mtx_t *W, FILE *f);
#endif
int avgw_write_info(avgw_mtx_t *W, FILE *f);
void avgw_hist2aw(grid_t *g, int n, int i0, int nsamp, const unsigned *h, avgw_func_t *f);
void avgw_itail(const grid_t *g, const avgw_func_t *f, const avgw_shapes_t *s, avgw_cutparam_t *c);
void avgw_reshape(const grid_t *g, const avgw_func_t *f, const gridparam_t *p, const avgw_cutparam_t *c, float *fcut);

#endif /* __RISM_AVGW_AVGW_H */
