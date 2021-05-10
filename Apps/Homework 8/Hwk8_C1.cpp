#include "Precompilied.h"
#include "Apps/HomeworkDrivers.h"

static constexpr int u = 0;

static real f(real x, real y) { return sinl(x + cosl(y)); }
static real dfdx(real x, real y) { return cosl(x + cosl(y)); }
static real dfdy(real x, real y) { return -1.0 * sinl(y) * cosl(x + cosl(y)); }

void Hwk8_C1_Driver()
{
  const int p = 2;
  const int n_gq = 7;
  const int K = 14;

  UnstructuredMesh2D mesh = UnstructuredMesh2D("Mesh.txt");
  FEM2D<1> fem = FEM2D<1>(mesh, p);
  L2_Projection2D(fem, f, n_gq);

  // Get coordinates of the center of element K
  const MeshNode2D& A1 = mesh(K, 0);
  const MeshNode2D& A2 = mesh(K, 1);
  const MeshNode2D& A3 = mesh(K, 2);
  const real x = (A1.x + A2.x + A3.x) / 3;
  const real y = (A1.y + A2.y + A3.y) / 3;

  // Calculate H1 error
  const real H1_x = FE_Error2DGlobal(u, fem, dfdx, n_gq, 1, 0);
  const real H1_y = FE_Error2DGlobal(u, fem, dfdy, n_gq, 0, 1);
  const real H1 = sqrtl(H1_x * H1_x + H1_y * H1_y);

  // Print values
  print(fem.evaluate(u, x, y));
  print(fem.evaluate(u, x, y, 1, 0));
  print(fem.evaluate(u, x, y, 0, 1));
  print(FE_Error2DGlobal(u, fem, f, n_gq, 0, 0));
  print(H1);
}
