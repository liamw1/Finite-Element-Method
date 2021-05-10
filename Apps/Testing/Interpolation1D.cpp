#pragma once
#include "Precompilied.h"
#include "Apps/HomeworkDrivers.h"

static real f(real x) { return expl(x) * sinl(2 * PI * x); }

void Interpolation_Driver()
{
  const real xL = 0, xR = 1;
  const int n = 3;

  for (int p = 1; p <= 3; ++p)
  {
    UniformMesh1D mesh = UniformMesh1D(xL, xR, n);
    FEM1D fem = FEM1D(mesh, p, f);

    std::string fileName = "Out";
    fileName += std::to_string(p);
    fileName += ".txt";
    fem.plot(fileName, false, 1000);
  }
}