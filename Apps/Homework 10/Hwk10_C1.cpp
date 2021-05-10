#include "Precompilied.h"
#include "Apps/HomeworkDrivers.h"

static constexpr int u1 = 0;
static constexpr int u2 = 1;
static constexpr int p = 0;

static real f1(real x, real y) { return x * (sinl(x) * sinl(y) + 3 * x) - y * cosl(x) * cosl(y) + 2 * (x * y + 1) * sinl(x) * cosl(y); }
static real f2(real x, real y) { return cosl(x) * (x * cosl(y) - 2 * (x * y + 1) * sinl(y)) - y * sinl(x) * sinl(y) - PI * cosl(PI * y); }
static real nu(real x, real y) { return 1 + x * y; }
static real rho(real x, real y) { return 1; }
static real g_D1(real x, real y) { return x * (x - 1) * (y - 1) * (cosl(y) - 1) + sinl(x) * cosl(y); }
static real g_D2(real x, real y) { return sinl(y) * ((x - 1) * (y - 1) * (cosl(x) - 1)) - cosl(x); }

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

void Hwk10_C1_Driver()
{
  const real xMin = 0, xMax = 1;
  const real yMin = 0, yMax = 1;
  const int nx = 80, ny = 80;
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

  // Print values
  const real x = PI / 7;
  const real y = PI / 5;
  print(uFem.evaluate(u1, x, y));
  print(uFem.evaluate(u2, x, y));
  print(pFem.evaluate(p, x, y));
  print();
  print(uFem.evaluate(u1, x, y, 1, 0));
  print(uFem.evaluate(u2, x, y, 1, 0));
  print(pFem.evaluate(p, x, y, 1, 0));
  print();
  print(uFem.evaluate(u1, x, y, 0, 1));
  print(uFem.evaluate(u2, x, y, 0, 1));
  print(pFem.evaluate(p, x, y, 0, 1));
  uFem.plot(u1);
  uFem.plot(u2);
  pFem.plot(p);
}