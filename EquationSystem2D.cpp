#include "Precompilied.h"
#include "EquationSystem2D.h"

Vector constructBoundaryVector2D(const FEM2D& fem, const Matrix& massMatrix)
{
  const int& p = fem.polynomialOrder;

  Vector bc_e = Vector(fem.Ng);
  for (int i = 0; i < fem.Ng; ++i)
    for (int n = 0; n < fem.boundaryIndices.size(); ++n)
    {
      const int& j = fem.boundaryIndices[n];
      bc_e[i] -= massMatrix[i][j] * fem.FENodes[j].u;
    }
  return bc_e;
}
