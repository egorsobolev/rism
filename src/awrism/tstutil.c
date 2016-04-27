#include <stdio.h>

void write_dmtx(const char *fn, int n, int *nn, double *d)
{
  FILE *f;
  int m, i;
  m = 1;
  for (i = 0; i < n; i++)
    m *= nn[i];
  f = fopen(fn, "wb");
  fwrite(&n, sizeof(int), 1, f);
  fwrite(nn, sizeof(int), n, f);
  fwrite(d, sizeof(double), m, f);
  fclose(f);
}

void write_imtx(const char *fn, int n, int *nn, int *d)
{
  FILE *f;
  int m, i;
  m = 1;
  for (i = 0; i < n; i++)
    m *= nn[i];
  f = fopen(fn, "wb");
  fwrite(&n, sizeof(int), 1, f);
  fwrite(nn, sizeof(int), n, f);
  fwrite(d, sizeof(int), m, f);
  fclose(f);
}
