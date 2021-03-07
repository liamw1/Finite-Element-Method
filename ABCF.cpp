#include "Precompilied.h"
#include "ABCF.h"

ABCF::ABCF(realFunction aFunc, realFunction bFunc, realFunction cFunc, realFunction fFunc, realFunction naturalBoundaryCondition)
  : a(aFunc), b(bFunc), c(cFunc), f(fFunc), naturalBC(naturalBoundaryCondition)
{
}

Vector ABCF::solveSystem(const FEM& fem, const int n_gq) const
{
  // Create mass matrices and load vector
  Matrix M_xx = FE_MassMatrix1D(fem, a, n_gq, 1);
  Matrix M_0x = FE_MassMatrix1D(fem, b, n_gq, 0, 1);
  Matrix M_00 = FE_MassMatrix1D(fem, c, n_gq, 0);
  Vector f_h = FE_LoadVector1D(fem, f, n_gq, 0);

  // Create boundary vectors
  Vector bc_n = constructNaturalBoundaryVector(fem, naturalBC);
  Vector bc_eM_xx = constructEssentialBoundaryVector(fem, M_xx);
  Vector bc_eM_0x = constructEssentialBoundaryVector(fem, M_0x);
  Vector bc_eM_00 = constructEssentialBoundaryVector(fem, M_00);

  // Remove boundary indices
  for (int n = (int)fem.boundaryIndices.size() - 1; n >= 0; --n)
  {
    const int& i = fem.boundaryIndices[n];
    M_xx.removeRowAndCol(i);
    M_0x.removeRowAndCol(i);
    M_00.removeRowAndCol(i);
    f_h.remove(i);
    bc_n.remove(i);
    bc_eM_xx.remove(i);
    bc_eM_0x.remove(i);
    bc_eM_00.remove(i);
  }

  // Solve linear system
  Vector coefficients = solve(M_xx + M_0x + M_00, f_h + bc_n + bc_eM_xx + bc_eM_0x + bc_eM_00);

  // Add back in boundary indices to coefficient vector
  for (int n = (int)fem.boundaryIndices.size() - 1; n >= 0; --n)
  {
    const int& i = fem.boundaryIndices[n];
    coefficients.insert(0, i);
  }

  return coefficients;
}
