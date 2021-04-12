#include "Precompilied.h"
#include "Integration.h"

real definiteIntegralSimpson(real1DFunction f, const real a, const real b)
{
  real h = (b - a) / N;
  real sum = 0.0;
  for (int i = 1; i <= N / 2; ++i)
    sum += f(a + 2.0 * (i - 1.0) * h) + 4 * f(a + (2.0 * i - 1.0) * h) + f(a + 2.0 * i * h);
  return sum * h / 3;
}