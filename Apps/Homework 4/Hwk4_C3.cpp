#include "Precompilied.h"
#include "Apps/HomeworkDrivers.h"

static real g_n(real x) { return (2.0 + expl(x)) * (cosl(x) - sinl(x)); }
static real a(real x) { return 2.0 + expl(x); }
static real b(real x) { return 0.0; }
static real c(real x) { return 1.0 + x * x; }
static real f(real x) { return (expl(x) * (x - 1.0) + x * (x * x + 3.0)) * cosl(x) + (4.0 + expl(x) * (x + 2.0)) * sinl(x); }
static real analyticalSolution(real x) { return x * cosl(x); }
static real analyticalSolutionDerivative(real x) { return cosl(x) - x * sinl(x); }

void Hwk4_C3_Driver()
{
  const real xL = 0, xR = 1;
  const int p = 2, n_gq = 3;

  std::ofstream file1("ErrorRegression.csv");
  std::ofstream file2("DerivativeErrorRegression.csv");
  for (int n = 10; n < 80; n += 10)
  {
    // Create mesh
    UniformMesh1D mesh = UniformMesh1D(xL, xR, n);
    mesh.setBoundaryConditions(BC_Type::Natural, BC_Type::Natural);
    FEM1D fem = FEM1D(mesh, p);

    // Create equation system
    Elliptic1DABCF equation = Elliptic1DABCF(fem, a, b, c, f, g_n);

    // Solve system
    equation.update(n_gq);
    print(FE_Error1DGlobal(fem, analyticalSolution, n_gq, 0));

    // Write regression
    const real h = (xR - xL) / n;
    file1 << logl(h) << ", " << logl(FE_Error1DGlobal(fem, analyticalSolution, n_gq, 0)) << std::endl;
    file2 << logl(h) << ", " << logl(FE_Error1DGlobal(fem, analyticalSolutionDerivative, n_gq, 1)) << std::endl;
  }
  file1.close();
  file2.close();
  std::cin.get();
  remove("ErrorRegression.csv");
  remove("DerivativeErrorRegression.csv");
}
