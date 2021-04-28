#include "Precompilied.h"
#include "HomeworkDrivers.h"

static constexpr int u = 0;

static real a(real x, real y) { return 1.0 + x + y / 2; }
static real b(real x, real y) { return 0.0; }
static real c(real x, real y) { return 2.0; }
static real g_N(real x, real y) { return 0.0; }
static real f(real x, real y)
{
  const real A = PI * cosl(2.0 * PI * y) * sinl(PI * x);
  const real B = 0.5 * cosl(PI * x);
  const real C = (4.0 + 5.0 * PI * PI * (2.0 + 2.0 * x + y)) * cosl(2.0 * PI * y);
  const real D = 2.0 * PI * sinl(2.0 * PI * y);
  return A + B * (C + D);
}

void Hwk9_C2_Driver()
{
  const real xMin = 0, xMax = 1;
  const real yMin = 0, yMax = 1;
  const int nx = 20, ny = 20;
  const int p = 1;
  const int n_gq = 7;

  // Create mesh and FEM structure
  UniformRectangularMesh2D mesh = UniformRectangularMesh2D(xMin, xMax, yMin, yMax, nx, ny);
  mesh.setBoundaryConditions(BC_Type::Natural);
  FEM2D<1> fem = FEM2D<1>(mesh, p);

  // Create equation system
  Elliptic2DABCF equationSystem = Elliptic2DABCF(fem, a, b, c, f, g_N);

  // Solve system
  equationSystem.update(n_gq);
  fem.plot(u, 1);

  // Print values
  const real x = PI / 4, y = PI / 6;
  print(fem.evaluate(u, x, y));
  print(fem.evaluate(u, x, y, 1, 0));
  print(fem.evaluate(u, x, y, 0, 1));
}