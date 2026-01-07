#include <stdlib.h>
#include <math.h>
#include <cblas.h>

typedef int Jx_func(void *eq_data, const float *x, float *r);

/* + 5 * n * sizeof(float)*/
int bicgstab(int N, const float *b, float *x, Jx_func *Jx, void *eq_data, float *tol, int *it, float *workspace)
{
	float *r, *h, *p, *v, *s;
	float rTh, rTr, normr, norms, alpha, beta, omega, sTr;
	float err, normb;
	int its, maxit;

	r = workspace;
	h = r + N;
	p = h + N;
	v = p + N;
	s = v + N;

	its = 0;
	maxit = *it * 2;
	normb = cblas_snrm2(N,b,1);
	if (normb == 0.0f)
		normb = 1.0f;
	err = *tol * normb;

	Jx(eq_data, x, r);
	++its;

	cblas_saxpy(N,-1.,b,1,r,1); /* r = Ax-b */
	cblas_scopy(N,r,1,h,1);
	cblas_scopy(N,r,1,p,1);

	rTh=cblas_sdot(N,r,1,h,1); /* rho1 */
	normr=cblas_snrm2(N,r,1);
	if (normr <= err) {
		*tol = normr / normb;
		*it = its;
		return 0;
	}
	while (its < maxit) {
		Jx(eq_data,p,v);  /* Jx(eq, p, v) */
		++its;

		alpha=rTh/cblas_sdot(N,h,1,v,1);
		cblas_scopy(N,r,1,s,1);  /* s = h */
		cblas_saxpy(N,-alpha,v,1,s,1); /* s = s - alpha * v */
		norms=cblas_snrm2(N,s,1);
		if (norms <= err) {
			cblas_saxpy(N,-alpha,p,1,x,1);
			*tol = norms / normb;
			*it = its;
			return 0;
		}

		Jx(eq_data,s,r);
		++its;
		sTr=cblas_sdot(N,s,1,r,1);
		rTr=cblas_sdot(N,r,1,r,1);
		if ( fabs(sTr)<1e-40 || fabs(rTr)<1e-40 )
			omega = 0.;
		else
			omega = sTr/rTr;

		cblas_saxpy(N,-alpha,p,1,x,1);
		cblas_saxpy(N,-omega,s,1,x,1);
		cblas_sscal(N,-omega,r,1);
		cblas_saxpy(N,1.f,s,1,r,1);

		normr=cblas_snrm2(N,r,1);
		if (normr <= err) {
			*tol = normr / normb;
			*it = its;
			return 0;
		}
		if (omega == 0.0f) {
			*tol = normr / normb;
			*it = its;
			return 2;
		}

		beta=(alpha/omega)/rTh;
		rTh=cblas_sdot(N,r,1,h,1);
		if (rTh == 0.0f) {
			*tol = normr / normb;
			*it = its;
			return 3;
		}
		beta*=rTh;

		cblas_sscal(N,beta,p,1);
		cblas_saxpy(N,1.,r,1,p,1);
		cblas_saxpy(N,-beta*omega,v,1,p,1);
	}

	*tol = rTr / normb;
	*it = its;
	return 1; /* no convergence */
}
