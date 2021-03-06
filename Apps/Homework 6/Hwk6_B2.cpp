#include "Precompilied.h"
#include "Apps/HomeworkDrivers.h"

static real g_n(real x) { return expl(x) * (3 * x + 1 + 2 * cosl(x)); }
static real a(real x) { return 2 + cosl(x); }
static real b(real x) { return 0; }
static real c(real x) { return 0; }
static real f(real x) { return expl(x) * ((x + 1) * sinl(x) - (x + 2) * (cosl(x) + 2)); }
static real analyticalSolution(real x) { return x * expl(x) - 0.7597602982176987; }
static real analyticalSoltionDerivative(real x) { return (1 + x) * expl(x); }

void Hwk6_B2_Driver()
{
  const real xL = 0, xR = 1;
  const int n = 20, p = 1, n_gq = 4;

  // Create mesh
  UniformMesh1D mesh = UniformMesh1D(xL, xR, n);
  mesh.setBoundaryConditions(BC_Type::Natural, BC_Type::Natural);
  FEM1D fem = FEM1D(mesh, p);

  // Create equation system
  Elliptic1DABCF equation = Elliptic1DABCF(fem, a, b, c, f, g_n);

  // Solve system
  equation.update(n_gq);
  print(fem.evaluate(PI / 6));
  print(FE_Error1DGlobal(fem, analyticalSolution, n_gq, 0));
}