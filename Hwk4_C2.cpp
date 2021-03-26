#include "Precompilied.h"
#include "Hwk4_C2.h"

static real g_e(real x) { return 1.0; }
static real g_n(real x) { return 9.0 * cosl(x) - 12.0 * sinl(x); }
static real a(real x) { return x + 2.0; }
static real b(real x) { return x + 1.0; }
static real c(real x) { return x * x + 1.0; }
static real f(real x) { return (3.0 * x * x * x + 4.0 * x * x + 13.0 * x + 3.0) * cosl(x) + (-3.0 * x * x + 5.0 * x + 12.0) * sinl(x); }

static void EnforceBoundaryConditions(FEM1D& fem)
{
  for (int K = 0; K < fem.meshSize; ++K)
    for (int j = 0; j < fem.polynomialOrder + 1; ++j)
      if (fem(K, j).BC == BC_Type::Essential)
        fem(K, j).u = g_e(fem(K, j).x);
}

void Hwk4_C2_Driver()
{
  const real xL = 0, xR = 1;
  const int n = 20, p = 2, n_gq = 3;

  // Create mesh
  UniformMesh1D mesh = UniformMesh1D(xL, xR, n);
  mesh.setBoundaryConditions(BC_Type::Essential, BC_Type::Natural);
  FEM1D fem = FEM1D(mesh, p);

  // Create equation system
  ABCF equation = ABCF(a, b, c, f, g_n);

  // Solve system
  EnforceBoundaryConditions(fem);
  update(fem, equation, n_gq);
  EnforceBoundaryConditions(fem);

  print(fem.evaluate(PI / 7));
  fem.plot(30, 0);
}
