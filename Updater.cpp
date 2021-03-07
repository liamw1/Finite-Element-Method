#include "Precompilied.h"
#include "Updater.h"

void update(FEM& fem, const EquationSystem1D& equationSystem, const int n_gq)
{
  Vector u_h = equationSystem.solveSystem(fem, n_gq);

  // Modify finite element solution
  for (int K = 0; K < fem.meshSize; ++K)
    for (int j = 0; j < fem.polynomialOrder + 1; ++j)
      fem(K, j).u = u_h[fem[K][j]];
}
