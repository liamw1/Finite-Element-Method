#include "Precompilied.h"
#include "LagrangeShapeFunctions2D.h"

real refLagrangePolynomial2D(const real tx, const real ty,
                             const int polynomialOrder, const int shapeIndex,
                             const int xDerivativeOrder, const int yDerivativeOrder)
{
  const int& j = shapeIndex;

  switch (polynomialOrder)
  {
  case 1:
    return refDegree1LagrangePolynomial2D(tx, ty, j, xDerivativeOrder, yDerivativeOrder);
  case 2:
    return refDegree2LagrangePolynomial2D(tx, ty, j, xDerivativeOrder, yDerivativeOrder);
  default:
    LOG("Polynomial order not supported", LogLevel::Error);
    return NAN;
  }
}

real refDegree1LagrangePolynomial2D(const real tx, const real ty,
                                    const int shapeIndex,
                                    const int xDerivativeOrder, const int yDerivativeOrder)
{
  ASSERT(xDerivativeOrder >= 0, "Derivative order in x must be non-negative");
  ASSERT(yDerivativeOrder >= 0, "Derivative order in y must be non-negative");
  ASSERT(shapeIndex >= 0, "Shape index must be non-negative");
  ASSERT(shapeIndex < 3, "Shape index must be less than 3");

  if (xDerivativeOrder == 0 && yDerivativeOrder == 0)
    switch (shapeIndex)
    {
    case 0:
      return 1 - tx - ty;
    case 1:
      return tx;
    case 2:
      return ty;
    }
  else if (xDerivativeOrder == 1 && yDerivativeOrder == 0)
    switch (shapeIndex)
    {
    case 0:
      return -1;
    case 1:
      return 1;
    case 2:
      return 0;
    }
  else if (xDerivativeOrder == 0 && yDerivativeOrder == 1)
    switch (shapeIndex)
    {
    case 0:
      return -1;
    case 1:
      return 0;
    case 2:
      return 1;
    }
  else
    return 0;

  LOG("No value was returned", LogLevel::Error);
  return NAN;
}

real refDegree2LagrangePolynomial2D(const real tx, const real ty,
                                    const int shapeIndex,
                                    const int xDerivativeOrder, const int yDerivativeOrder)
{
  ASSERT(xDerivativeOrder >= 0, "Derivative order in x must be non-negative");
  ASSERT(yDerivativeOrder >= 0, "Derivative order in y must be non-negative");
  ASSERT(shapeIndex >= 0, "Shape index must be non-negative");
  ASSERT(shapeIndex < 6, "Shape index must be less than 6");

  const real phi0 = refDegree1LagrangePolynomial2D(tx, ty, 0, 0, 0);
  const real phi1 = refDegree1LagrangePolynomial2D(tx, ty, 1, 0, 0);
  const real phi2 = refDegree1LagrangePolynomial2D(tx, ty, 2, 0, 0);

  if (xDerivativeOrder == 0 && yDerivativeOrder == 0)
    switch (shapeIndex)
    {
    case 0:
      return 2 * phi0 * (phi0 - 0.5);
    case 1:
      return 2 * phi1 * (phi1 - 0.5);
    case 2:
      return 2 * phi2 * (phi2 - 0.5);
    case 3:
      return 4 * phi0 * phi1;
    case 4:
      return 4 * phi2 * phi1;
    case 5:
      return 4 * phi2 * phi0;
    }
  else if (xDerivativeOrder == 1 && yDerivativeOrder == 0)
  {
    const real phi0x = refDegree1LagrangePolynomial2D(tx, ty, 0, 1, 0);
    const real phi1x = refDegree1LagrangePolynomial2D(tx, ty, 1, 1, 0);
    const real phi2x = refDegree1LagrangePolynomial2D(tx, ty, 2, 1, 0);

    switch (shapeIndex)
    {
    case 0:
      return phi0x * (4 * phi0 - 1);
    case 1:
      return phi1x * (4 * phi1 - 1);
    case 2:
      return phi2x * (4 * phi2 - 1);
    case 3:
      return 4 * (phi0x * phi1 + phi0 * phi1x);
    case 4:
      return 4 * (phi2x * phi1 + phi2 * phi1x);
    case 5:
      return 4 * (phi2x * phi0 + phi2 * phi0x);
    }
  }
  else if (xDerivativeOrder == 0 && yDerivativeOrder == 1)
  {
    const real phi0y = refDegree1LagrangePolynomial2D(tx, ty, 0, 0, 1);
    const real phi1y = refDegree1LagrangePolynomial2D(tx, ty, 1, 0, 1);
    const real phi2y = refDegree1LagrangePolynomial2D(tx, ty, 2, 0, 1);

    switch (shapeIndex)
    {
    case 0:
      return phi0y * (4 * phi0 - 1);
    case 1:
      return phi1y * (4 * phi1 - 1);
    case 2:
      return phi2y * (4 * phi2 - 1);
    case 3:
      return 4 * (phi0y * phi1 + phi0 * phi1y);
    case 4:
      return 4 * (phi2y * phi1 + phi2 * phi1y);
    case 5:
      return 4 * (phi2y * phi0 + phi2 * phi0y);
    }
  }
  else if (xDerivativeOrder == 1 && yDerivativeOrder == 1)
  {
    const real phi0x = refDegree1LagrangePolynomial2D(tx, ty, 0, 1, 0);
    const real phi1x = refDegree1LagrangePolynomial2D(tx, ty, 1, 1, 0);
    const real phi2x = refDegree1LagrangePolynomial2D(tx, ty, 2, 1, 0);
    const real phi0y = refDegree1LagrangePolynomial2D(tx, ty, 0, 0, 1);
    const real phi1y = refDegree1LagrangePolynomial2D(tx, ty, 1, 0, 1);
    const real phi2y = refDegree1LagrangePolynomial2D(tx, ty, 2, 0, 1);

    switch (shapeIndex)
    {
    case 0:
      return 4 * phi0x * phi0y;
    case 1:
      return 4 * phi1x * phi1y;
    case 2:
      return 4 * phi2x * phi2y;
    case 3:
      return 4 * (phi0x * phi1y + phi0y * phi1x);
    case 4:
      return 4 * (phi1x * phi2y + phi1y * phi2x);
    case 5:
      return 4 * (phi0x * phi2y + phi0y * phi2x);
    }
  }
  else if (xDerivativeOrder == 2 && yDerivativeOrder == 0)
  {
    const real phi0x = refDegree1LagrangePolynomial2D(tx, ty, 0, 1, 0);
    const real phi1x = refDegree1LagrangePolynomial2D(tx, ty, 1, 1, 0);
    const real phi2x = refDegree1LagrangePolynomial2D(tx, ty, 2, 1, 0);

    switch (shapeIndex)
    {
    case 0:
      return 4 * phi0x * phi0x;
    case 1:
      return 4 * phi1x * phi1x;
    case 2:
      return 4 * phi2x * phi2x;
    case 3:
      return 8 * phi0x * phi1x;
    case 4:
      return 8 * phi1x * phi2x;
    case 5:
      return 8 * phi0x * phi2x;
    }
  }
  else if (xDerivativeOrder == 0 && yDerivativeOrder == 2)
  {
    const real phi0y = refDegree1LagrangePolynomial2D(tx, ty, 0, 0, 1);
    const real phi1y = refDegree1LagrangePolynomial2D(tx, ty, 1, 0, 1);
    const real phi2y = refDegree1LagrangePolynomial2D(tx, ty, 2, 0, 1);

    switch (shapeIndex)
    {
    case 0:
      return 4 * phi0y * phi0y;
    case 1:
      return 4 * phi1y * phi1y;
    case 2:
      return 4 * phi2y * phi2y;
    case 3:
      return 8 * phi0y * phi1y;
    case 4:
      return 8 * phi1y * phi2y;
    case 5:
      return 8 * phi0y * phi2y;
    }
  }
  else
    return 0;

  LOG("No value was returned", LogLevel::Error);
  return NAN;
}

void plotRefLagrangePolynomial(const int degree, const int shapeIndex, const int xDerivativeOrder, const int yDerivativeOder, const int n)
{
  std::ofstream file("Plot.txt");
  for (int i = 0; i <= n; ++i)
    for (int j = 0; j <= n - i; ++j)
    {
      const real x = (real)i / n;
      const real y = (real)j / n;
      real f = 0.0;

      switch (degree)
      {
      case 1:
        f = refDegree1LagrangePolynomial2D(x, y, shapeIndex, xDerivativeOrder, yDerivativeOder);
        break;
      case 2:
        f = refDegree2LagrangePolynomial2D(x, y, shapeIndex, xDerivativeOrder, yDerivativeOder);
        break;
      default:
        LOG("Polynomial degree not implemented", LogLevel::Error);
      }

      file << x << ", " << y << ", " << f << std::endl;
    }
  file.close();
  std::cin.get();
  remove("Plot.txt");
}

template <int N>
void plotShapeFunction(const FEM2D<N>& fem, const int elementIndex, const int shapeIndex, const int xDerivativeOrder, const int yDerivativeOrder, const int n)
{
  const int& K = elementIndex;
  const MeshNode2D& A1 = fem.mesh(K, 0);
  const MeshNode2D& A2 = fem.mesh(K, 1);
  const MeshNode2D& A3 = fem.mesh(K, 2);

  // Create transformation matrix
  Matrix B = Matrix(2);
  B[0][0] = A2.x - A1.x;
  B[0][1] = A3.x - A1.x;
  B[1][0] = A2.y - A1.y;
  B[1][1] = A3.y - A1.y;

  // Plot polynomial
  std::ofstream file("Plot.txt");
  for (int i = 0; i <= n; ++i)
    for (int j = 0; j <= n - i; ++j)
    {
      // Create reference points
      const real tx = (real)i / n;
      const real ty = (real)j / n;

      // Transform points from reference domain to element K
      const real x = B[0][0] * tx + B[0][1] * ty + A1.x;
      const real y = B[1][0] * tx + B[1][1] * ty + A1.y;

      real f = lagrangeShapeFunction2D(x, y, fem, K, shapeIndex, xDerivativeOrder, yDerivativeOrder);
      file << x << ", " << y << ", " << f << std::endl;
    }
  file.close();
  std::cout << "Data written to \"Plot.txt\".  Press any key to continue" << std::endl;
  std::cin.get();
  remove("Plot.txt");
}
