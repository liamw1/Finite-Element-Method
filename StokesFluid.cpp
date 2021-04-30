#include "Precompilied.h"
#include "StokesFluid.h"

static constexpr int u1 = 0;
static constexpr int u2 = 1;
static constexpr int p = 0;
static real identityFunction(real x, real y) { return 1.0; };

StokesFluid::StokesFluid(FEM2D<2>& uFEM, FEM2D<1>& pFEM,
                         real2DFunction f1Func, real2DFunction f2Func,
                         real2DFunction nuFunc,
                         real2DFunction rhoFunc)
  : uFem(uFEM), pFem(pFEM),
    f1(f1Func), f2(f2Func),
    nu(nuFunc),
    rho(rhoFunc)
{
}

int StokesFluid::neq() const
{
  return 3;
}

Vector StokesFluid::solveSystem(const int n_gq) const
{
  std::function<real(real, real)> rhoInverse = [=](real x, real y) { return 1.0 / rho(x, y); };

  // Create mass matrices and load vector
  Matrix Muu_xx = FE_MassMatrix2D(uFem, nu, n_gq, 1, 0);
  Matrix Muu_yy = FE_MassMatrix2D(uFem, nu, n_gq, 0, 1);
  Matrix Mpu_0x = FE_MassMatrix2D(pFem, uFem, rhoInverse, n_gq, 0, 0, 1, 0);
  Matrix Mpu_0y = FE_MassMatrix2D(pFem, uFem, rhoInverse, n_gq, 0, 0, 0, 1);
  Vector f_l = FE_LoadVector2D(pFem, identityFunction, n_gq, 0, 0);
  Vector f1u_h = FE_LoadVector2D(uFem, f1, n_gq, 0, 0);
  Vector f2u_h = FE_LoadVector2D(uFem, f2, n_gq, 0, 0);

  // Create boundary vectors
  Vector bc_u1_eM_xx = constructEssentialBoundaryVector2D(uFem, u1, Muu_xx);
  Vector bc_u1_eM_yy = constructEssentialBoundaryVector2D(uFem, u1, Muu_yy);
  Vector bc_u1_eM_0x = constructEssentialBoundaryVector2D(pFem, uFem, u1, Mpu_0x);

  Vector bc_u2_eM_xx = constructEssentialBoundaryVector2D(uFem, u2, Muu_xx);
  Vector bc_u2_eM_yy = constructEssentialBoundaryVector2D(uFem, u2, Muu_yy);
  Vector bc_u2_eM_0y = constructEssentialBoundaryVector2D(pFem, uFem, u2, Mpu_0y);

  // Remove boundary indices (since we already know the values at those nodes)
  removeBoundaryIndices(Muu_xx, uFem.boundaryIndices);
  removeBoundaryIndices(Muu_yy, uFem.boundaryIndices);
  removeBoundaryIndices(Mpu_0x, pFem.boundaryIndices, uFem.boundaryIndices);
  removeBoundaryIndices(Mpu_0y, pFem.boundaryIndices, uFem.boundaryIndices);
  removeBoundaryIndices(f_l, pFem.boundaryIndices);
  removeBoundaryIndices(f1u_h, uFem.boundaryIndices);
  removeBoundaryIndices(f2u_h, uFem.boundaryIndices);
  removeBoundaryIndices(bc_u1_eM_xx, uFem.boundaryIndices);
  removeBoundaryIndices(bc_u1_eM_yy, uFem.boundaryIndices);
  removeBoundaryIndices(bc_u1_eM_0x, pFem.boundaryIndices);
  removeBoundaryIndices(bc_u2_eM_xx, uFem.boundaryIndices);
  removeBoundaryIndices(bc_u2_eM_yy, uFem.boundaryIndices);
  removeBoundaryIndices(bc_u2_eM_0y, pFem.boundaryIndices);

  // Construct full coefficient matrix and boundary vector
  const int Nu_u = uFem.Ng - (int)uFem.boundaryIndices.size();  // Number of unknowns on each component of u
  const int Nu_p = pFem.Ng - (int)pFem.boundaryIndices.size();  // Number of unknowns on p

  // Debug
  ASSERT(Muu_xx.size() == Nu_u, "Matrix is not the correct size");
  ASSERT(Muu_yy.size() == Nu_u, "Matrix is not the correct size");
  ASSERT(Mpu_0x.rows() == Nu_p && Mpu_0x.columns() == Nu_u, "Matrix is not the correct size");
  ASSERT(Mpu_0y.rows() == Nu_p && Mpu_0x.columns() == Nu_u, "Matrix is not the correct size");
  ASSERT(f_l.size() == Nu_p, "Vector is not the correct size");
  ASSERT(f1u_h.size() == Nu_u, "Vector is not the correct size");
  ASSERT(f2u_h.size() == Nu_u, "Vector is not the correct size");
  ASSERT(bc_u1_eM_xx.size() == Nu_u, "Vector is not the correct size");
  ASSERT(bc_u1_eM_yy.size() == Nu_u, "Vector is not the correct size");
  ASSERT(bc_u1_eM_0x.size() == Nu_p, "Vector is not the correct size");
  ASSERT(bc_u2_eM_xx.size() == Nu_u, "Vector is not the correct size");
  ASSERT(bc_u2_eM_yy.size() == Nu_u, "Vector is not the correct size");
  ASSERT(bc_u2_eM_0y.size() == Nu_p, "Vector is not the correct size");

  // Matrix M = Matrix(2 * Nu_u + Nu_p + 1);
  std::vector<std::vector<real>> M = std::vector<std::vector<real>>();
  M.resize(2 * Nu_u + Nu_p + 1);
  for (int i = 0; i < M.size(); ++i)
    M[i].resize(2 * Nu_u + Nu_p + 1, 0.0);

  for (int i = 0; i < Nu_u; ++i)
    for (int j = 0; j < Nu_u; ++j)
    {
      M[i][j] = Muu_xx[i][j] + Muu_yy[i][j];
      M[Nu_u + i][Nu_u + j] = Muu_xx[i][j] + Muu_yy[i][j];
    }
  for (int i = 0; i < Nu_p; ++i)
    for (int j = 0; j < Nu_u; ++j)
    {
      M[2 * Nu_u + i][j] = -1.0 * Mpu_0x[i][j];
      M[2 * Nu_u + i][Nu_u + j] = -1.0 * Mpu_0y[i][j];
      M[j][2 * Nu_u + i] = -1.0 * Mpu_0x[i][j];
      M[Nu_u + j][2 * Nu_u + i] = -1.0 * Mpu_0y[i][j];
    }
  for (int i = 0; i < Nu_p; ++i)
  {
    M[2 * Nu_u + Nu_p][2 * Nu_u + i] = f_l[i];
    M[2 * Nu_u + i][2 * Nu_u + Nu_p] = f_l[i];
  }

  Vector b = Vector(2 * Nu_u + Nu_p + 1);
  for (int i = 0; i < Nu_u; ++i)
  {
    b[i] = f1u_h[i] + bc_u1_eM_xx[i] + bc_u1_eM_yy[i];
    b[Nu_u + i] = f2u_h[i] + bc_u2_eM_xx[i] + bc_u2_eM_yy[i];
  }
  for (int i = 0; i < Nu_p; ++i)
    b[2 * Nu_u + i] = -1.0 * (bc_u1_eM_0x[i] + bc_u2_eM_0y[i]);

  // Solve linear system for coefficients on unknown nodes
  Vector coefficients = solve(M, b);

  // Add back in boundary indices to coefficient vector
  for (int n = 0; n < uFem.boundaryIndices.size(); ++n)
  {
    const int& i = uFem.boundaryIndices[n];
    coefficients.insert(0, i);
  }
  for (int n = 0; n < uFem.boundaryIndices.size(); ++n)
  {
    const int& i = uFem.boundaryIndices[n];
    coefficients.insert(0, uFem.Ng + i);
  }
  for (int n = 0; n < pFem.boundaryIndices.size(); ++n)
  {
    const int& i = pFem.boundaryIndices[n];
    coefficients.insert(0, 2 * uFem.Ng + i);
  }

  // Debug
  ASSERT(coefficients.size() == 2 * uFem.Ng + pFem.Ng + 1, "Coefficient vector is not the correct size");
 
  return coefficients;
}

void StokesFluid::update(const int n_gq)
{
  Vector u_h = solveSystem(n_gq);

  // Modify finite element solution
  for (int i = 0; i < uFem.Ng; ++i)
  {
    uFem.FENodes[i][u1] = u_h[i];
    uFem.FENodes[i][u2] = u_h[uFem.Ng + i];
  }
  for (int i = 0; i < pFem.Ng; ++i)
    pFem.FENodes[i][p] = u_h[2 * uFem.Ng + i];
}
