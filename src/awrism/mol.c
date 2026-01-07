#include "mol.h"

#define _CRT_SECURE_NO_WARNINGS

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int mol_read(const char *fn, mol_t *m)
{
	FILE *f;
	char b[20];
	int n, c, i, j, l;
	static int fl[] = {12, 12, 12, 10, 16, 10};
	double *p[6];
	int *uniq;

	f = fopen(fn, "rt");
	if (!f) {
		return 1;
	}
	n = 0;
	while (!feof(f) && !n) {
		fgets(b, 9, f);
		l = strlen(b);
		if (b[0] != '#' && l != 0) {
			n = atoi(b);
			if (n < 1)
				return 2;
		}
		if (l == 8) {
			c = fgetc(f);
			/* skip rest of line */
			while (!feof(f) && c != '\n')
				c = fgetc(f);
		}
	}
	m->x = p[0] = (double *) malloc(n * (3 * sizeof(double) + sizeof(int)));
	if (!m->x)
		return 3;
	m->y = p[1] = m->x + n;
	m->z = p[2] = m->y + n;
	m->atyp = (int *) (m->z + n);
	p[3] = (double *) malloc(n * (3 * sizeof(double) + sizeof(int)));
	if (!p[3]) {
		free(m->x);
		return 3;
	}
	p[4] = p[3] + n;
	p[5] = p[4] + n;
	uniq = (int *) (p[5] + n);

	i = 0;
	while (!feof(f) && i < n) {
		c = fgetc(f);
		if (c != '#' && c != '\n') {
			ungetc(c, f);
			for (j = 0; j < 6; ++j) {
				fgets(b, fl[j] + 1, f);
				l = strlen(b);
				if (l < fl[j]) {
					free(m->x);
					free(p[3]);
					return 2;
				}
				p[j][i] = atof(b);
			}
			c = (char) fgetc(f);
			++i;
		}
		/* skip rest of line */
		while (!feof(f) && c != '\n')
			c = (char) fgetc(f);
	}
	fclose(f);
	if (i < n) {
		free(m->x);
		free(p[3]);
		return 2;
	}

	m->natom = n;

	m->qtot = 0;
	for (i = 0; i < n; i++)
		m->qtot += p[4][i];

	m->ntype = mol_uniq(m->natom, p[3], uniq, m->atyp);

	m->ff.s = (double *) malloc(m->ntype * 3 * sizeof(double));
	if (!m->ff.s) {
		free(m->x);
		free(p[3]);
		return 4;
	}
	m->ff.q = m->ff.s + m->ntype;
	m->ff.e = m->ff.q + m->ntype;

	for (i = 0; i < m->ntype; i++) {
		j = uniq[i];
		m->ff.s[i] = p[3][j];
		m->ff.q[i] = p[4][j];
		m->ff.e[i] = p[5][j];
	}
	free(p[3]);

	return 0;
}

void mol_del(mol_t *m)
{
	free(m->x);
	free(m->ff.s);
}

int mol_uniq(int n, const double *ff, int *uniq, int *orig)
{
	int i, j, k, l, notexist, nn;

	nn = 2 * n;
	uniq[0] = 0;
	orig[0] = 0;
	j = 1;
	for (i = 1; i < n; i++) {
		notexist = (ff[i] != ff[0]) || (ff[i+n] != ff[n]) || (ff[i+nn] != ff[nn]);
		k = 1;
		while (k < j && notexist) {
			l = uniq[k];
			notexist = (ff[i] != ff[l]) || (ff[i+n] != ff[l+n]) || (ff[i+nn] != ff[l+nn]);
			k++;
		}
		if (notexist) {
			uniq[j] = i;
			j++;
		} else
			k--;

		orig[i] = k;
	}
	return j;
}
