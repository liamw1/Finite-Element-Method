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
  uFem.plotField();
  pFem.plot(p);
  
  // Print values
  const real x1 = -0.9, y1 = 0.9;
  const real x2 = 0.9, y2 = 0.9;
  print(pFem.evaluate(p, x1, y1));
  print(pFem.evaluate(p, x2, y2));
  print();

  // Use Newton's method to find u(x, y) = 0
  const int maxIterations = 100;
  real x = 0.0, y = 0.5;
  for (int n = 0; n < maxIterations; ++n)
  {
    // Compute Jacobian
    Matrix J = Matrix(2);
    J[0][0] = uFem.evaluate(u1, x, y, 1, 0);  J[0][1] = uFem.evaluate(u1, x, y, 0, 1);
    J[1][0] = uFem.evaluate(u2, x, y, 1, 0);  J[1][1] = uFem.evaluate(u2, x, y, 0, 1);

    // Invert matrix
    const real determinant = J[0][0] * J[1][1] - J[0][1] * J[1][0];
    ASSERT(determinant != 0.0, "Transformation matrix is singular");
    real temp = J[0][0];
    J[0][0] = J[1][1];
    J[1][1] = temp;
    J[0][1] *= -1;
    J[1][0] *= -1;
    J /= determinant;

    // Update solution
    const real eps1 = J[0][0] * uFem.evaluate(u1, x, y) + J[0][1] * uFem.evaluate(u2, x, y);
    const real eps2 = J[1][0] * uFem.evaluate(u1, x, y) + J[1][1] * uFem.evaluate(u2, x, y);
    x -= eps1;
    y -= eps2;

    if (x < xMin || x > xMax || y < yMin || y > yMax)
      LOG("(x, y) is not in the domain", LogLevel::Error);
  }
  print("(", x, ", ", y, ")");
  print("(", uFem.evaluate(u1, x, y), ", ", uFem.evaluate(u2, x, y), ")");
}