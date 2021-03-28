#include "Hwk4_C1.h"
#include "Hwk4_C2.h"
#include "Hwk4_C3.h"
#include "Hwk6_B2.h"
#include "UniformRectangularMesh2D.h"
#include "LagrangeShapeFunctions2D.h"

int main()
{
  const int K = 2;
  const int j = 5;
  const int d_x = 0;
  const int d_y = 0;
  const int n = 20;

  UniformRectangularMesh2D mesh = UniformRectangularMesh2D(0, 1, 0.5, 1.1, 2, 2);
  FEM2D fem = FEM2D(mesh, 2);

  print("(", mesh(K, 0).x, ", ", mesh(K, 0).y, "), (", mesh(K, 1).x, ", ", mesh(K, 1).y, "), (", mesh(K, 2).x, ", ", mesh(K, 2).y, ")");
  plotShapeFunction(fem, K, j, d_x, d_y, n);
}