#include <math.h>
#include <float.h>

int fractostr(char *s, double a, double eps)
{
  double r;
  int d, i;
  if (eps == 0.0)
    eps = DBL_EPSILON;
  else
    eps = fabs(eps);
  r = fabs(a - (int) (a + eps));
  i = 0;
  while (r > eps) {
    r *= 10;
    eps *= 10;
    d = (int) (r + eps);
    s[i] = d ^ 0x30;
    r -= d;
    i++;
  }
  return i;
}
