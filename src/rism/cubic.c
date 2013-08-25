#include <math.h>

void sspline(int n1, const float *x, const float *y, float yp1, float ypn, float *y2, float *u)
{
  int i, n;
  float sig, p, qn, un;

  n = n1 - 1;

  if (isinf(yp1)) {
    y2[0] = u[0] = 0.0; 
  } else {
    y2[0] = -0.5;
    u[0] = (3.0 / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1);
  }

  for (i = 1; i < n; i++) {
    sig = (x[i] - x[i-1]) / (x[i+1] - x[i-1]);
    p = sig * y2[i-1] + 2.0;
    y2[i] = (sig - 1.0) / p;
    u[i] = (y[i+1] - y[i]) / (x[i+1] - x[i]) - (y[i] - y[i-1]) / (x[i] - x[i-1]);
    u[i] = (6.0 * u[i] / (x[i+1] - x[i-1]) - sig * u[i-1]) / p;
  }
  if (isinf(ypn)) {
    qn = un = 0.0;
  } else {
    qn = 0.5;
    un = (3.0 / (x[n] - x[n-1])) * (ypn - (y[n] - y[n-1]) / (x[n] - x[n-1]));
  }

  y2[n] = (un - qn * u[n-1]) / (qn * y2[n-1] + 1.0);
  for (i = n-1; i >= 0; i--)
    y2[i] = y2[i] * y2[i+1] + u[i]; 
}

void ssplint(int n, const float *xa, const float *ya, const float *y2a, int m, float *xx, float *yy)
{
  int klo, khi, i, k;
  float h, b, a, x;

  for (i = 0; i < m; i++) {
    klo = 1;
    khi = n;

    x = xx[i];
    while ((khi - klo) > 1) {
      k = (khi + klo) >> 1;
      if (xa[k] > x) khi = k;
      else klo = k;
    }
    h = xa[khi] - xa[klo];
    if (h == 0.0) {
      yy[i] = NAN;
    } else {
      a = (xa[khi] - x) / h; 
      b = (x - xa[klo]) / h; 
      yy[i] = a * ya[klo] + b * ya[khi] + ((a*a*a-a) * y2a[klo] + (b*b*b-b) * y2a[khi]) * (h*h) / 6.0;
    }
  }
}

void sspline_uni(int n1, const float *y, float yp1, float ypn, float *y2, float *u)
{
  int i, n;
  float sig, p, qn, un;

  n = n1 - 1;

  if (isinf(yp1)) {
    y2[0] = u[0] = 0.0; 
  } else {
    y2[0] = -0.5;
    u[0] = 3.0 * (y[1] - y[0] - yp1);
  }

  for (i = 1; i < n; i++) {
    y2[i] = -1.0 / (y2[i-1] + 4.0);
    u[i] = (y[i+1] - y[i]) - (y[i] - y[i-1]);
    u[i] = (6.0 * u[i] - u[i-1]) / (y2[i-1] + 4.0);
  }
  if (isinf(ypn)) {
    qn = un = 0.0;
  } else {
    qn = 0.5;
    un = 3.0 * (ypn - y[n] + y[n-1]);
  }

  y2[n] = (un - qn * u[n-1]) / (qn * y2[n-1] + 1.0);
  for (i = n-1; i >= 0; i--)
    y2[i] = y2[i] * y2[i+1] + u[i]; 
}

void ssplint_uni(int n, const float *ya, const float *y2a, int m, float x0, float h, float *yy)
{
  float x, a, b;
  int klo, khi, i;

  x = x0;
  for (i = 0; i < m; i++) {
    klo = (int) floor(x);
    khi = (int) ceil(x);
    if (klo < 0) {
      klo = 0;
      khi = 1;
    } else if (khi >= n) {
      khi = n - 1;
      klo = n - 2;
    }

    a = khi - x;
    b = x - klo;
    yy[i] = a * ya[klo] + b * ya[khi] + ((a*a*a-a) * y2a[klo] + (b*b*b-b) * y2a[khi]) / 6.0;
 
    x += h;
  }
}
