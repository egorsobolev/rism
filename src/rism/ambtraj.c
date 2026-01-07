#include <rism/ambtraj.h>

#include <errno.h>
#include <stdio.h>

int amb_open(const char *fname, int natm, int ewald, ambtraj_t *traj)
{
	int c, i, n;
	traj->f = fopen(fname, "r");
	if (!traj->f)
		return -1;
	i = 0;
	c = fgetc(traj->f);
	while (c != EOF && c != '\n' && i < AMBTRAJ_TITLE_LEN) {
		traj->title[i] = c;
		i++;
		c = fgetc(traj->f);
	}
	traj->title[i] = 0;
	while (c != EOF && c != '\n')
		c = fgetc(traj->f);

	i = 3 * natm;
	n = i / 10;
	if (i - n * 10)
		n++;
	if (ewald)
		n++;

	traj->err = 0;
	traj->frm = 0;
	traj->natm = natm;
	traj->ewald = ewald;
	traj->nline = n;

	return 0;
}
int amb_read(ambtraj_t *traj, float *x)
{
	int n, i;
	float a, b, c;
	if (traj->err)
		return -1;

	n = 3 * traj->natm;

	if (fscanf(traj->f, "%f", x) != 1)
		return 0;
	for (i = 1; i < n; i++) {
		if (fscanf(traj->f, "%f", x + i) != 1) {
			traj->err = errno;
			return -1;
		}
	}
	if (traj->ewald) {
		if (fscanf(traj->f, "%f%f%f", &a, &b, &c) != 3) {
			traj->err = errno;
			return -1;
		}
	}
	i = fgetc(traj->f);
	while (i != EOF && i != '\n')
		i = fgetc(traj->f);
	traj->frm++;
	return 1;
}
int amb_ignore(ambtraj_t *traj, int nskip)
{
	int i, c, n;
	if (traj->err)
		return -1;

	n = 0;
	while (!feof(traj->f) && n < nskip) {
		i = 0;
		c = fgetc(traj->f);
		if (c != EOF) {
			do {
				if (c == '\n')
					i++;
			} while (i < traj->nline && (c = fgetc(traj->f)) != EOF);
			if (i < traj->nline) {
				traj->err = errno;
				return -1;
			}
			traj->frm++;
			n++;
		}
	}
	return n;
}
void amb_close(ambtraj_t *traj)
{
	fclose(traj->f);
}
