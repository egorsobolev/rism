#ifndef __RISM_GRID_H
#define __RISM_GRID_H

#include <stdlib.h>

struct Grid
{
  size_t n, ngrid;
  double dr, dk, f, b, c;
};
typedef struct Grid grid_t;

int grid_init(int n, double dr, grid_t *g);

#endif /* __RISM_GRID_H */
