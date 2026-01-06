#ifndef __RISM_GRID_H
#define __RISM_GRID_H

struct GRID_FLOAT
{
  int np;
  float dr;
  float dt;
  float f;
  float b;
};

typedef struct GRID_FLOAT sgrid_t;

struct GRID_DOUBLE
{
  int np;
  double dr;
  double dt;
  double f;
  double b;
};

typedef struct GRID_DOUBLE dgrid_t;

int sgrid_init(int np, float dr, sgrid_t *g);
int dgrid_init(int np, double dr, dgrid_t *g);

#endif /* __RISM_GRID_H */
