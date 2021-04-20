#include "Precompilied.h"
#include "HomeworkDrivers.h"

static real a(real x, real y) { return 1.0; }
static real b(real x, real y) { return x / 9; }
static real c(real x, real y) { return 2.0; }
static real f(real x, real y) { return (-18.0 + x * x * (20.0 + 3 * y * y) + x * y * (-54.0 + 19.0 * y * y)) / 9; }
static real g_D(real x, real y) { return x * x + y * y * y * (x - sinl(PI * x)) + y * sinl(PI * x); }

void EnforceBoundaryConditions(FEM2D& fem)
{
  for (int i = 0; i < fem.Ng; ++i)
    if (fem.FENodes[i].BC == BC_Type::Dirichlet || fem.FENodes[i].BC == BC_Type::Corner)
    {
      const real& x = fem.FENodes[i].x;
      const real& y = fem.FENodes[i].y;
      fem.FENodes[i].u = g_D(x, y);
    }
}

void Hwk9_C1_Driver()
{
  const real xMin = 0, xMax = 1;
  const real yMin = 0, yMax = 1;
  const int nx = 40, ny = 40;
  const int p = 1;
  const int n_gq = 7;

  UniformRectangularMesh2D mesh = UniformRectangularMesh2D(xMin, xMax, yMin, yMax, nx, ny);
  mesh.setBoundaryConditions(BC_Type::Dirichlet);
  FEM2D fem = FEM2D(mesh, p);
  Elliptic2DABCF equationSystem = Elliptic2DABCF(a, b, c, f);

  EnforceBoundaryConditions(fem);
  update2D(fem, equationSystem, n_gq);
  EnforceBoundaryConditions(fem);

  const int x = PI / 4, y = PI / 6;
  fem.plot(1);
}