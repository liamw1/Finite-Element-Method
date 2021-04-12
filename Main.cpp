#include "Hwk4_C1.h"
#include "Hwk4_C2.h"
#include "Hwk4_C3.h"
#include "Hwk6_B2.h"
#include "L2Projection.h"
#include "UniformRectangularMesh2D.h"

real f(real x, real y)
{
  return cosl(x + sinl(y));
}

real dfdx(real x, real y)
{
  return -1.0 * sinl(x + sinl(y));
}

real dfdy(real x, real y)
{
  return -1.0 * cosl(y) * sinl(x + sinl(y));
}


int main()
{
  const int p = 2;
  const int n_gq = 7;

  UniformRectangularMesh2D mesh1 = UniformRectangularMesh2D(0, 1, 0, 1, 40, 40);
  FEM2D fem1 = FEM2D(mesh1, p, f);

  const real err_L2_1 = FE_Error2DGlobal(fem1, f, 7, 0, 0);
  const real err_H1_x_1 = FE_Error2DGlobal(fem1, dfdx, 7, 1, 0);
  const real err_H1_y_1 = FE_Error2DGlobal(fem1, dfdy, 7, 0, 1);
  const real err_H1_1 = sqrtl(err_H1_x_1 * err_H1_x_1 + err_H1_y_1 * err_H1_y_1);


  UniformRectangularMesh2D mesh2 = UniformRectangularMesh2D(0, 1, 0, 1, 80, 80);
  FEM2D fem2 = FEM2D(mesh2, p, f);

  const real err_L2_2 = FE_Error2DGlobal(fem2, f, 7, 0, 0);
  const real err_H1_x_2 = FE_Error2DGlobal(fem2, dfdx, 7, 1, 0);
  const real err_H1_y_2 = FE_Error2DGlobal(fem2, dfdy, 7, 0, 1);
  const real err_H1_2 = sqrtl(err_H1_x_2 * err_H1_x_2 + err_H1_y_2 * err_H1_y_2);


  print(err_L2_1 / err_L2_2);
  print(err_H1_1 / err_H1_2);
}