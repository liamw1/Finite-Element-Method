#pragma once
#include "Precompilied.h"
#include "Meshing/2D/FEM2D.h"
#include "LinearAlgebra/Matrix.h"

/*
  Calculates a p-th degree Lagrange polynomial on the
  reference triangle [(0, 0), (0, 1), (1, 0)].

  \param tx: Result when x is mapped onto reference domain
  \param ty: Result when y is mapped onto reference domain
*/
real refLagrangePolynomial2D(const real tx, const real ty,
                             const int polynomialOrder, const int shapeIndex,
                             const int xDerivativeOrder, const int yDerivativeOrder);

/*
  Calculates a first degree Lagrange polynomial on the
  reference triangle [(0, 0), (0, 1), (1, 0)].

  \param tx: Result when x is mapped onto reference domain
  \param ty: Result when y is mapped onto reference domain
*/
real refDegree1LagrangePolynomial2D(const real tx, const real ty, 
                                    const int shapeIndex,
                                    const int xDerivativeOrder, const int yDerivativeOrder);

/*
  Calculates a first degree Lagrange polynomial on the
  reference triangle [(0, 0), (0, 1), (1, 0)].

  \param tx: Result when x is mapped onto reference domain
  \param ty: Result when y is mapped onto reference domain
*/
real refDegree2LagrangePolynomial2D(const real tx, const real ty,
                                    const int shapeIndex,
                                    const int xDerivativeOrder, const int yDerivativeOrder);

void plotRefLagrangePolynomial(const int degree, const int shapeIndex, const int xDerivativeOrder, const int yDerivativeOder, const int n);

template<int N>
void plotShapeFunction(const FEM2D<N>& fem,
                       const int elementIndex, const int shapeIndex,
                       const int xDerivativeOrder, const int yDerivativeOrder,
                       const int n);

template<int N>
real lagrangeShapeFunction2D(const real x, const real y,
  const FEM2D<N>& fem,
  const int elementIndex, const int shapeIndex,
  const int xDerivativeOrder, const int yDerivativeOrder)
{
  const int& K = elementIndex;
  const int& j = shapeIndex;
  const MeshNode2D& A1 = fem.mesh(K, 0);
  const MeshNode2D& A2 = fem.mesh(K, 1);
  const MeshNode2D& A3 = fem.mesh(K, 2);

  // Create transformation matrix from reference domain to K
  Matrix B = Matrix(2);
  B[0][0] = A2.x - A1.x;  B[0][1] = A3.x - A1.x;
  B[1][0] = A2.y - A1.y;  B[1][1] = A3.y - A1.y;

  // Invert matrix
  const real determinant = B[0][0] * B[1][1] - B[0][1] * B[1][0];
  ASSERT(determinant != 0.0, "Transformation matrix is singular");
  real temp = B[0][0];
  B[0][0] = B[1][1];
  B[1][1] = temp;
  B[0][1] *= -1;
  B[1][0] *= -1;
  B /= determinant;

  const real tx = B[0][0] * (x - A1.x) + B[0][1] * (y - A1.y);
  const real ty = B[1][0] * (x - A1.x) + B[1][1] * (y - A1.y);
  ASSERT(tx > -1.0 * TOLERANCE && ty > -1.0 * TOLERANCE && ty < 1.0 - tx + TOLERANCE, "(tx, ty) is not inside the reference domain");

  if (xDerivativeOrder == 0 && yDerivativeOrder == 0)
    return refLagrangePolynomial2D(tx, ty, fem.polynomialOrder, j, 0, 0);
  else if (xDerivativeOrder == 1 && yDerivativeOrder == 0)
  {
    return B[0][0] * refLagrangePolynomial2D(tx, ty, fem.polynomialOrder, j, 1, 0)
      + B[1][0] * refLagrangePolynomial2D(tx, ty, fem.polynomialOrder, j, 0, 1);
  }
  else if (xDerivativeOrder == 0 && yDerivativeOrder == 1)
  {
    return B[0][1] * refLagrangePolynomial2D(tx, ty, fem.polynomialOrder, j, 1, 0)
      + B[1][1] * refLagrangePolynomial2D(tx, ty, fem.polynomialOrder, j, 0, 1);
  }
  else
  {
    LOG("Derivative order not implemented", LogLevel::Error);
    return 0.0;
  }
}