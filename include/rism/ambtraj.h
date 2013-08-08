#ifndef __RISM_AMBTRAJ_H
#define __RISM_AMBTRAJ_H

#include <stdio.h>

#define AMBTRAJ_TITLE_LEN 80

struct AMBTRAJ
{
  FILE *f;
  int natm;
  int ewald;
  int err;
  int frm;
  int nline;
  char title[AMBTRAJ_TITLE_LEN + 1];
};

typedef struct AMBTRAJ ambtraj_t;

int amb_open(const char *fname, int n, int ewald, ambtraj_t *traj);
int amb_read(ambtraj_t *traj, float *x);
int amb_ignore(ambtraj_t* traj, int nskip);
void amb_close(ambtraj_t *traj);


#endif /* __RISM_AMBTRAJ_H */
