#pragma once
#include "Precompilied.h"
#include "EquationSystem1D.h"
#include "EquationSystem2D.h"

/*
  Updates solution based on equation system provided.

  Note: Essential boundary conditions should be enforeced
  before each call of this function.
*/
void update1D(FEM1D& fem, const EquationSystem1D& equationSystem, const int n_gq);

/*
  Updates solution based on equation system provided.

  Note: Essential boundary conditions should be enforeced
  before each call of this function.
*/
template<int N>
void update2D(FEM2D<N>& fem, const EquationSystem2D& equationSystem, const int n_gq)
{
  const int& p = fem.polynomialOrder;
  Vector u_h = equationSystem.solveSystem(n_gq);

  // Modify finite element solution
  for (int K = 0; K < fem.mesh.size; ++K)
    for (int j = 0; j < (p + 1) * (p + 2) / 2; ++j)
    {
      ASSERT(!std::isinf(u_h[fem[K][j]]), "FE update results in Infinite value");
      ASSERT(!std::isnan(u_h[fem[K][j]]), "FE update results in NaN");
      fem(K, j)[0] = u_h[fem[K][j]];
    }
}