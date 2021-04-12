#include "Precompilied.h"
#include "EquationSytem1D.h"

Vector constructNaturalBoundaryVector(const FEM1D& fem, real1DFunction naturalBC)
{
  Vector bc_n = Vector(fem.Ng);
  for (int K = 0; K < fem.meshSize; ++K)
    for (int j = 0; j < fem.polynomialOrder + 1; ++j)
    {
      if (fem(K, j).BC == BC_Type::Natural)
      {
        bc_n[fem[K][j]] = naturalBC(fem(K, j).x);
        if (fem(K, j).x == fem.mesh.xL)  // Invert sign for left natural boundary
          bc_n[fem[K][j]] *= -1.0;
      }

      // Skip last node in K to avoid double-counting
      if (K != fem.meshSize - 1 && j == fem.polynomialOrder - 1)
        ++j;
    }
  return bc_n;
}

Vector constructEssentialBoundaryVector(const FEM1D& fem, const Matrix& massMatrix)
{
  Vector bc_e = Vector(fem.Ng);
  for (int K = 0; K < fem.meshSize; ++K)
    for (int j = 0; j < fem.polynomialOrder + 1; ++j)
    {
      // Using essential boundary conditions
      for (int n = 0; n < fem.boundaryIndices.size(); ++n)
      {
        const int& i = fem.boundaryIndices[n];
        bc_e[fem[K][j]] -= massMatrix[fem[K][j]][i] * fem.FENodes[i].u;
      }

      // Skip last node in K to avoid double-counting
      if (K != fem.meshSize - 1 && j == fem.polynomialOrder - 1)
        ++j;
    }
  return bc_e;
}
