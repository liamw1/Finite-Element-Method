#include "Precompilied.h"
#include "Updater.h"

void update1D(FEM1D& fem, const EquationSystem1D& equationSystem, const int n_gq)
{
  Vector u_h = equationSystem.solveSystem(fem, n_gq);

  // Modify finite element solution
  for (int K = 0; K < fem.meshSize; ++K)
    for (int j = 0; j < fem.polynomialOrder + 1; ++j)
    {
      ASSERT(!std::isinf(u_h[fem[K][j]]), "FE update results in Infinite value");
      ASSERT(!std::isnan(u_h[fem[K][j]]), "FE update results in NaN");
      fem(K, j).u = u_h[fem[K][j]];
    }
}

void update2D(FEM2D& fem, const EquationSystem2D& equationSystem, const int n_gq)
{
  const int& p = fem.polynomialOrder;
  Vector u_h = equationSystem.solveSystem(fem, n_gq);

  // Modify finite element solution
  for (int K = 0; K < fem.mesh.size; ++K)
    for (int j = 0; j < (p + 1) * (p + 2) / 2; ++j)
    {
      ASSERT(!std::isinf(u_h[fem[K][j]]), "FE update results in Infinite value");
      ASSERT(!std::isnan(u_h[fem[K][j]]), "FE update results in NaN");
      fem(K, j).u = u_h[fem[K][j]];
    }
}
