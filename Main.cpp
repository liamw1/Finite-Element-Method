#include "Hwk4_C1.h"
#include "Hwk4_C2.h"
#include "Hwk4_C3.h"
#include "Hwk6_B2.h"
#include "UniformRectangularMesh2D.h"
#include "LagrangeShapeFunctions2D.h"

int main()
{
  const int p = 3;
  const int nx = 1;
  const int ny = 1;

  UniformRectangularMesh2D mesh = UniformRectangularMesh2D(0, p * nx , 0, p * ny, nx, ny);
  FEM2D fem = FEM2D(mesh, p);

  for (int i = 0; i < fem.Ng; ++i)
    print("(", fem.FENodes[i].x, ", ", fem.FENodes[i].y, ")");
}