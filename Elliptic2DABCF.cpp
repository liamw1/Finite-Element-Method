#include "Precompilied.h"
#include "Elliptic2DABCF.h"

static constexpr int u = 0;

Elliptic2DABCF::Elliptic2DABCF(FEM2D<1>& uFem, real2DFunction aFunc, real2DFunction bFunc, real2DFunction cFunc, real2DFunction fFunc, real2DFunction naturalBoundaryCondition)
  : fem(uFem), a(aFunc), b(bFunc), c(cFunc), f(fFunc), naturalBC(naturalBoundaryCondition)
{
}

int Elliptic2DABCF::neq() const
{
  return 1;
}

Vector Elliptic2DABCF::solveSystem(const int n_gq) const
{
  // Create mass matrices and load vector
  Matrix* M_xx = FE_MassMatrix2D(fem, a, n_gq, 1, 0);
  Matrix* M_yy = FE_MassMatrix2D(fem, a, n_gq, 0, 1);
  Matrix* M_0x = FE_MassMatrix2D(fem, b, n_gq, 0, 0, 1, 0);
  Matrix* M_0y = FE_MassMatrix2D(fem, b, n_gq, 0, 0, 0, 1);
  Matrix* M_00 = FE_MassMatrix2D(fem, c, n_gq, 0, 0);
  Vector* f_h = FE_LoadVector2D(fem, f, n_gq, 0, 0);

  // Create boundary vectors
  Vector* bc_n = constructNaturalBoundaryVector2D(fem, naturalBC, n_gq);
  Vector* bc_eM_xx = constructEssentialBoundaryVector2D(fem, u, M_xx);
  Vector* bc_eM_yy = constructEssentialBoundaryVector2D(fem, u, M_yy);
  Vector* bc_eM_0x = constructEssentialBoundaryVector2D(fem, u, M_0x);
  Vector* bc_eM_0y = constructEssentialBoundaryVector2D(fem, u, M_0y);
  Vector* bc_eM_00 = constructEssentialBoundaryVector2D(fem, u, M_00);

  // Remove boundary indices (since we already know the values at those nodes)
  removeBoundaryIndices(M_xx, fem.boundaryIndices);
  removeBoundaryIndices(M_yy, fem.boundaryIndices);
  removeBoundaryIndices(M_0x, fem.boundaryIndices);
  removeBoundaryIndices(M_0y, fem.boundaryIndices);
  removeBoundaryIndices(M_00, fem.boundaryIndices);
  removeBoundaryIndices(f_h, fem.boundaryIndices);
  removeBoundaryIndices(bc_n, fem.boundaryIndices);
  removeBoundaryIndices(bc_eM_xx, fem.boundaryIndices);
  removeBoundaryIndices(bc_eM_yy, fem.boundaryIndices);
  removeBoundaryIndices(bc_eM_0x, fem.boundaryIndices);
  removeBoundaryIndices(bc_eM_0y, fem.boundaryIndices);
  removeBoundaryIndices(bc_eM_00, fem.boundaryIndices);

  // Solve linear system for coefficients on unknown nodes
  Vector coefficients = solve(*M_xx + *M_yy + *M_0x + *M_0y + *M_00, *f_h + *bc_n + *bc_eM_xx + *bc_eM_yy + *bc_eM_0x + *bc_eM_0y + *bc_eM_00);

  // Add back in boundary indices to coefficient vector
  for (int n = 0; n < fem.boundaryIndices.size(); ++n)
  {
    const int& i = fem.boundaryIndices[n];
    coefficients.insert(0, i);
  }

  // Free memory
  delete M_xx;
  delete M_yy;
  delete M_0x;
  delete M_0y;
  delete M_00;
  delete f_h;
  delete bc_n;
  delete bc_eM_xx;
  delete bc_eM_yy;
  delete bc_eM_0x;
  delete bc_eM_0y;
  delete bc_eM_00;

  return coefficients;
}

void Elliptic2DABCF::update(const int n_gq)
{
  const int& p = fem.polynomialOrder;
  Vector u_h = solveSystem(n_gq);

  // Modify finite element solution
  for (int K = 0; K < fem.mesh.size; ++K)
    for (int j = 0; j < (p + 1) * (p + 2) / 2; ++j)
    {
      ASSERT(!std::isinf(u_h[fem[K][j]]), "FE update results in Infinite value");
      ASSERT(!std::isnan(u_h[fem[K][j]]), "FE update results in NaN");
      fem(K, j)[u] = u_h[fem[K][j]];
    }
}
