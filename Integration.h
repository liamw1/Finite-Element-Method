#pragma once
#include "Precompilied.h"

const int N = 100;

/*
  Calculates definite integral of f from a to b using Simpson's rule.
*/
real definiteIntegralSimpson(real1DFunction f, const real a, const real b);