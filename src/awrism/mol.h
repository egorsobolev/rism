#ifndef __RISM_MOL_H
#define __RISM_MOL_H

struct MolForceFiled {
  double *s, *q, *e;
} ff;
typedef struct MolForceFiled molff_t;


struct Molecule {
  int natom, ntype;
  int *atyp, *ityp;
  double *x, *y, *z;
  molff_t ff;
  double o[3];
  double qtot;
};
typedef struct Molecule mol_t;


int mol_read(const char *, mol_t *);
int mol_uniq(int, const double *, int *, int *);
void mol_del(mol_t *);

#endif //__RISM_MOL_H
