#include "Precompilied.h"
#include "ACF.h"

ACF::ACF(real1DFunction aFunc, real1DFunction cFunc, real1DFunction fFunc, real1DFunction naturalBoundaryCondition)
  : a(aFunc), c(cFunc), f(fFunc), naturalBC(naturalBoundaryCondition)
{
}

Vector ACF::solveSystem(const FEM1D& fem, const int n_gq) const
{
  // Create mass matrices and load vector
  Matrix M_xx = FE_MassMatrix1D(fem, a, n_gq, 1);
  Matrix M_00 = FE_MassMatrix1D(fem, c, n_gq, 0);
  Vector f_h = FE_LoadVector1D(fem, f, n_gq, 0);

  // Create boundary vectors
  Vector bc_n = constructNaturalBoundaryVector(fem, naturalBC);
  Vector bc_eM_xx = constructEssentialBoundaryVector(fem, M_xx);
  Vector bc_eM_00 = constructEssentialBoundaryVector(fem, M_00);

  // Remove boundary indices
  for (int n = (int)fem.boundaryIndices.size() - 1; n >= 0; --n)
  {
    const int& i = fem.boundaryIndices[n];
    M_xx.removeRowAndCol(i);
    M_00.removeRowAndCol(i);
    f_h.remove(i);
    bc_n.remove(i);
    bc_eM_xx.remove(i);
    bc_eM_00.remove(i);
  }

  // Solve linear system
  Vector coefficients = solve(M_xx + M_00, f_h + bc_n + bc_eM_xx + bc_eM_00);

  // Add back in boundary indices to coefficient vector
  for (int n = (int)fem.boundaryIndices.size() - 1; n >= 0; --n)
  {
    const int& i = fem.boundaryIndices[n];
    coefficients.insert(0, i);
  }

  return coefficients;
}
