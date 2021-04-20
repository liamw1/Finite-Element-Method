#include "Precompilied.h"
#include "Elliptic2DABCF.h"

Elliptic2DABCF::Elliptic2DABCF(real2DFunction aFunc, real2DFunction bFunc, real2DFunction cFunc, real2DFunction fFunc)
  : a(aFunc), b(bFunc), c(cFunc), f(fFunc)
{
}

Vector Elliptic2DABCF::solveSystem(const FEM2D& fem, const int n_gq) const
{
  // Create mass matrices and load vector
  Matrix M_xx = FE_MassMatrix2D(fem, a, n_gq, 1, 0);
  Matrix M_yy = FE_MassMatrix2D(fem, a, n_gq, 0, 1);
  Matrix M_0x = FE_MassMatrix2D(fem, b, n_gq, 0, 0, 1, 0);
  Matrix M_0y = FE_MassMatrix2D(fem, b, n_gq, 0, 0, 0, 1);
  Matrix M_00 = FE_MassMatrix2D(fem, c, n_gq, 0, 0);
  Vector f_h = FE_LoadVector2D(fem, f, n_gq, 0, 0);

  // Create boundary vectors
  Vector bc_eM_xx = constructBoundaryVector2D(fem, M_xx);
  Vector bc_eM_yy = constructBoundaryVector2D(fem, M_yy);
  Vector bc_eM_0x = constructBoundaryVector2D(fem, M_0x);
  Vector bc_eM_0y = constructBoundaryVector2D(fem, M_0y);
  Vector bc_eM_00 = constructBoundaryVector2D(fem, M_00);

  // Remove boundary indices
  for (int n = (int)fem.boundaryIndices.size() - 1; n >= 0; --n)
  {
    const int& i = fem.boundaryIndices[n];
    M_xx.removeRowAndCol(i);
    M_yy.removeRowAndCol(i);
    M_0x.removeRowAndCol(i);
    M_0y.removeRowAndCol(i);
    M_00.removeRowAndCol(i);
    f_h.remove(i);
    bc_eM_xx.remove(i);
    bc_eM_yy.remove(i);
    bc_eM_0x.remove(i);
    bc_eM_0y.remove(i);
    bc_eM_00.remove(i);
  }

  // Solve linear system
  Vector coefficients = solve(M_xx + M_yy + M_0x + M_0y + M_00, f_h + bc_eM_xx + bc_eM_yy + bc_eM_0x + bc_eM_0y + bc_eM_00);

  // Add back in boundary indices to coefficient vector
  for (int n = 0; n < fem.boundaryIndices.size(); ++n)
  {
    const int& i = fem.boundaryIndices[n];
    coefficients.insert(0, i);
  }

  return coefficients;
}
