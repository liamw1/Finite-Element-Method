#pragma once
#include "Precompilied.h"
#include "FEM2D.h"
#include "Matrix.h"

real lagrangeShapeFunction2D(const real x, const real y,
                             const FEM2D& fem,
                             const int elementIndex, const int shapeIndex,
                             const int xDerivativeOrder, const int yDerivativeOrder);

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

void plotShapeFunction(const FEM2D& fem,
                       const int elementIndex, const int shapeIndex,
                       const int xDerivativeOrder, const int yDerivativeOrder,
                       const int n);