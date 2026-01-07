#ifndef __RISM_NR_H
#define __RISM_NR_H


typedef int eq_func(void *eq_data, const double *x, double *z, double *en);
typedef int Jx_func(void *eq_data, const float *x, float *r);
typedef void getx_func(const void *eq_data, const double *z, float *x);
typedef void putx_func(const void *eq_data, const float *x, double *z);

struct NRPARM
{
	int maxi, lmaxi, narm;
	double rtol, etol;
};
typedef struct NRPARM nrparm_t;

struct EQUATION
{
	int nZ;
	int nJx;
	eq_func *Z;
	Jx_func *Jx;
	getx_func *getb;
	putx_func *putx;
	void *p;
	nrparm_t nrprm;
	float *nr_solver_data;
	float *linear_solver_data;
};
typedef struct EQUATION eq_t;

int nr(eq_t *eq, double *t, double *rtol, double *etol, int *maxit);
int bicgstab(int N, const float *b, float *x, Jx_func *Jx, void *eq_data, float *tol, int *it, float *workspace);

#endif //__RISM_NR_H
