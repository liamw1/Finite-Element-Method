#include "Precompilied.h"
#include "HomeworkDrivers.h"

static constexpr int u1 = 0;
static constexpr int u2 = 1;
static constexpr int p = 0;

static real f1(real x, real y) { return 0; }
static real f2(real x, real y) { return 0; }
static real nu(real x, real y) { return 1; }
static real rho(real x, real y) { return 1; }
static real g_D1(real x, real y) { return (y > 1 - TOLERANCE && y < 1 + TOLERANCE) ? 1 : 0; }
static real g_D2(real x, real y) { return 0; }

static void EnforceBoundaryConditions(FEM2D<2>& fem)
{
  for (int i = 0; i < fem.Ng; ++i)
    if (fem.FENodes[i].BC == BC_Type::Dirichlet)
    {
      const real& x = fem.FENodes[i].x;
      const real& y = fem.FENodes[i].y;
      fem.FENodes[i][0] = g_D1(x, y);
      fem.FENodes[i][1] = g_D2(x, y);
    }
}

void Hwk10_C2_Driver()
{
  const real xMin = -1, xMax = 1;
  const real yMin = -1, yMax = 1;
  const int nx = 20, ny = 20;
  const int polyOrder = 1;
  const int n_gq = 7;

  // Create mesh and FEM structure
  UniformRectangularMesh2D mesh1 = UniformRectangularMesh2D(xMin, xMax, yMin, yMax, nx, ny);
  UniformRectangularMesh2D mesh2 = UniformRectangularMesh2D(xMin, xMax, yMin, yMax, nx, ny);
  mesh1.setBoundaryConditions(BC_Type::Dirichlet);
  mesh2.setBoundaryConditions(BC_Type::Natural);
  FEM2D<2> uFem = FEM2D<2>(mesh1, polyOrder + 1);
  FEM2D<1> pFem = FEM2D<1>(mesh2, polyOrder);

  // Create equation system
  StokesFluid eq = StokesFluid(uFem, pFem, f1, f2, nu, rho);

  // Solve system
  EnforceBoundaryConditions(uFem);
  eq.update(n_gq);
  EnforceBoundaryConditions(uFem);

  // Plot fields
  uFem.plot(u1);
  uFem.plot(u2);
  pFem.plot(p);
}