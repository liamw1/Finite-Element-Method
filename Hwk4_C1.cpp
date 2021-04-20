#include "Precompilied.h"
#include "Hwk4_C1.h"

static real g_e(real x) { return expl(x); }
static real g_n(real x) { return 0.0; }
static real a(real x) { return x * x + 2.0; }
static real c(real x) { return x + 1.0; }
static real f(real x) { return (-4.0 * x * x * x * x - 14.0 * x * x + x - 3.0) *   expl(x * x); }

static void EnforceBoundaryConditions(FEM1D& fem)
{
  for (int K = 0; K < fem.meshSize; ++K)
    for (int j = 0; j < fem.polynomialOrder + 1; ++j)
      if (fem(K, j).BC == BC_Type::Essential)
        fem(K, j).u = g_e(fem(K, j).x);
}

void Hwk4_C1_Driver()
{
  const real xL = 0, xR = 1;
  const int n = 20, p = 2, n_gq = 3;

  // Create mesh
  UniformMesh1D mesh = UniformMesh1D(xL, xR, n);
  mesh.setBoundaryConditions(BC_Type::Natural, BC_Type::Essential);
  FEM1D fem = FEM1D(mesh, p);

  // Create equation system
  Elliptic1DACF equation = Elliptic1DACF(a, c, f, g_n);

  // Solve system
  EnforceBoundaryConditions(fem);
  update1D(fem, equation, n_gq);
  EnforceBoundaryConditions(fem);

  print(fem.evaluate(PI / 4, 1));
  print(fem.integrate(n_gq));
  fem.plot(30, 0);
}
