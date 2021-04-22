#include "Precompilied.h"
#include "Elliptic2DACF.h"

Elliptic2DACF::Elliptic2DACF(real2DFunction aFunc, real2DFunction cFunc, real2DFunction fFunc, real2DFunction naturalBoundaryCondition)
  : a(aFunc), c(cFunc), f(fFunc), naturalBC(naturalBoundaryCondition)
{
}

Vector Elliptic2DACF::solveSystem(const FEM2D& fem, const int n_gq) const
{
  // Create mass matrices and load vector
  Matrix M_xx = FE_MassMatrix2D(fem, a, n_gq, 1, 0);
  Matrix M_yy = FE_MassMatrix2D(fem, a, n_gq, 0, 1);
  Matrix M_00 = FE_MassMatrix2D(fem, c, n_gq, 0, 0);
  Vector f_h = FE_LoadVector2D(fem, f, n_gq, 0, 0);

  // Create boundary vectors
  Vector bc_n = constructNaturalBoundaryVector2D(fem, naturalBC, n_gq);
  Vector bc_eM_xx = constructEssentialBoundaryVector2D(fem, M_xx);
  Vector bc_eM_yy = constructEssentialBoundaryVector2D(fem, M_yy);
  Vector bc_eM_00 = constructEssentialBoundaryVector2D(fem, M_00);

  // Remove boundary indices
  removeBoundaryIndices(M_xx, fem.boundaryIndices);
  removeBoundaryIndices(M_yy, fem.boundaryIndices);
  removeBoundaryIndices(M_00, fem.boundaryIndices);
  removeBoundaryIndices(f_h, fem.boundaryIndices);
  removeBoundaryIndices(bc_n, fem.boundaryIndices);
  removeBoundaryIndices(bc_eM_xx, fem.boundaryIndices);
  removeBoundaryIndices(bc_eM_yy, fem.boundaryIndices);
  removeBoundaryIndices(bc_eM_00, fem.boundaryIndices);

  // Solve linear system
  Vector coefficients = solve(M_xx + M_yy + M_00, f_h + bc_n + bc_eM_xx + bc_eM_yy + bc_eM_00);

  // Add back in boundary indices to coefficient vector
  for (int n = 0; n < fem.boundaryIndices.size(); ++n)
  {
    const int& i = fem.boundaryIndices[n];
    coefficients.insert(0, i);
  }

  return coefficients;
}
