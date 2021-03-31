#include "Hwk4_C1.h"
#include "Hwk4_C2.h"
#include "Hwk4_C3.h"
#include "Hwk6_B2.h"
#include "UniformRectangularMesh2D.h"
#include "LagrangeShapeFunctions2D.h"

real f(real x, real y)
{
  return sinl(x) + cosl(y);
}


int main()
{
  const int p = 1;
  const int nx = 100;
  const int ny = 100;
  const int n_gq = 7;

  UniformRectangularMesh2D mesh = UniformRectangularMesh2D(0, 1 , 0, 1, nx, ny);
  FEM2D fem = FEM2D(mesh, p);

  const int K = 0;

  real sum = 0.0;
  for (int K = 0; K < mesh.size; ++K)
  {
    const auto t = gauss2DNodesLocal(mesh, K, n_gq);
    const auto w = gauss2DWeightsLocal(mesh, K, n_gq);
    for (int i = 0; i < n_gq; ++i)
      sum += w[i] * f(t[i][0], t[i][1]);
  }
  print(sum);
}