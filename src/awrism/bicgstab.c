#include <stdlib.h>
#include <math.h>

typedef int Jx_func(void *eq_data, const float *x, float *r);

int bicgstab(int N, const float *b, float *x, Jx_func *Jx, void *eq_data, float *tol, int *it, float *workspace)
{
	float *r, *h, *p, *v, *s;
	float rTh, rTr, sTr, hTv;
	float sumb2, sumr2, sums2;
	int its, maxit, more, exitcode;

	/*
	 r, h, p, v, s: float[N]
	 */
	r = workspace;
	h = r + N;
	p = h + N;
	v = p + N;
	s = v + N;

	sumb2 = 0.0f;
	rTh = 0.0f;
	sumr2 = 0.0f;
	#pragma omp parallel
	{
		int i, maxit, its, more;
		float normb, normr, norms;
		float err, omega, beta, alpha;

		maxit = *it * 2;
		its = 0;
		more = 1;

		#pragma omp for reduction(+:sumb2)
		for (i = 0; i < N; i++)
			sumb2 += b[i] * b[i];
		normb = sqrtf(sumb2);
		if (normb == 0.0f)
			normb = 1.0f;
		err = *tol * normb;
		Jx(eq_data, x, r);

		++its;

		#pragma omp for reduction(+:rTh, sumr2)
		for (i = 0; i < N; i++) {
			r[i] -= b[i];
			p[i] = h[i] = r[i];
			rTh += r[i] * h[i];
			sumr2 += r[i] * r[i];
		}
		normr = sqrtf(sumr2);
		if (normr <= err) {
			#pragma omp single
			{
				*tol = normr / normb;
				*it = its;
				exitcode = 0;
			}
			more = 0;
		}
		while (more && its < maxit) {
			Jx(eq_data, p, v);  /* Jx(eq, p, v) */
			++its;
			#pragma omp single nowait
			hTv = 0.0f;
			#pragma omp for reduction(+:hTv)
			for (i = 0; i < N; i++)
				hTv += h[i] * v[i];
			alpha = rTh / hTv;
			#pragma omp single nowait
			sums2 = 0.0f;
			#pragma omp for reduction(+:sums2)
			for (i = 0; i < N; i++) {
				s[i] = r[i] - alpha * v[i];
				sums2 += s[i];
			}
			norms = sqrtf(sums2);
			if (norms <= err) {
				#pragma omp for
				for (i = 0; i < N; i++)
					x[i] -= alpha * p[i];
				#pragma omp single
				{
					*tol = norms / normb;
					*it = its;
					exitcode = 0;
				}
				more = 0;
				break;
			}
			Jx(eq_data, s, r);
			++its;
			#pragma omp single nowait
			sTr = 0.0f;
			#pragma omp single nowait
			rTr = 0.0f;
			#pragma omp for reduction(+:sTr, rTr)
			for (i = 0; i < N; i++) {
				sTr += s[i] * r[i];
				rTr += r[i] * r[i];
			}
			if ( fabs(sTr)<1e-40 || fabs(rTr)<1e-40 )
				omega = 0.;
			else
				omega = sTr/rTr;
			#pragma omp single nowait
			sumr2 = 0.0f;
			#pragma omp for reduction(+:sumr2)
			for (i = 0; i < N; i++) {
				x[i] -= alpha * p[i] + omega * s[i];
				r[i] = s[i] -omega * r[i];
				sumr2 += r[i] * r[i];
			}
			normr = sqrtf(sumr2);
			if (normr <= err) {
				#pragma omp single
				{
					*tol = normr / normb;
					*it = its;
					exitcode = 0;
				}
				more = 0;
				break;
			}
			if (omega == 0.0f) {
				#pragma omp single
				{
					*tol = normr / normb;
					*it = its;
					exitcode = 2;
				}
				more = 0;
				break;
			}
			beta=(alpha/omega)/rTh;
			#pragma omp barrier
			#pragma omp single nowait
			rTh = 0.0f;
			#pragma omp for reduction(+:rTh)
			for (i = 0; i < N; i++)
				rTh += r[i] * h[i];
			if (rTh == 0.0f) {
				#pragma omp single
				{
					*tol = normr / normb;
					*it = its;
					exitcode = 3;
				}
				more = 0;
				break;
			}
			beta*=rTh;
			#pragma omp for
			for (i = 0; i < N; i++)
				p[i] = beta * p[i] + r[i] - beta * omega * v[i];
		}
		if (more) {
			#pragma omp single
			{
				*tol = rTr / normb;
				*it = its;
				exitcode = 1;
			}
		}
	}
	return exitcode; /* no convergence */
}
