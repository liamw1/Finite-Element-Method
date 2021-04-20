#pragma once
#include "Precompilied.h"
#include "Matrix.h"
#include "Gauss-LegendreNodes.h"
#include "LagrangeShapeFunctions1D.h"
#include "LagrangeShapeFunctions2D.h"

/*
  \returns the FE load vector for a function f.

  \param f: Function in the inner products of the load vector.
  \param n_gq: Number of Gaussian quadrature nodes.
  \param derivativeOrder: Order of derivative on the Lagrange shape functions.  
  
  Note: This function does NOT cacluate the derivative of f.  
  Must pass in the derivative manually.
*/
Vector FE_LoadVector1D(const FEM1D& fem, real1DFunction f, const int n_gq, const int derivativeOrder);

/*
  \returns the FE mass matrix for a function "a" using one FEM1D.

  \param a: Function in the inner products of the mass matrix.
  \param n_gq: Number of Gaussian quadrature nodes.
  \param derivativeOrder: Order of derivative on the Lagrange shape functions.

  Note: Full matrix representation is inefficient here.
  Better to use sparse matrix representation.
*/
Matrix FE_MassMatrix1D(const FEM1D& fem, real1DFunction a, const int n_gq, const int derivativeOrder);

/*
  \returns the FE mass matrix for a function "a" using one FEM1D.

  \param a: Function in the inner products of the mass matrix.
  \param n_gq: Number of Gaussian quadrature nodes.
  \param derivativeOrder1: Order of derivative on the Lagrange shape functions on test functions v.
  \param derivativeOrder2: Order of derivative on the Lagrange shape functions on trial functions u.

  Note: Full matrix representation is inefficient here.
  Better to use sparse matrix representation.
*/
Matrix FE_MassMatrix1D(const FEM1D& fem, real1DFunction a, const int n_gq, const int derivativeOrder1, const int derivativeOrder2);

/*
  Performs and L2 projection on FEM1D for a function f.

  \param n_gq: Number of Gaussian quadrature nodes.
*/
void L2_Projection1D(FEM1D& fem, real1DFunction f, const int n_gq);



/*
  \returns the FE load vector for afunction f.

  \param f: Function in the inner products of the load vector.
  \param n_gq: Number of Gaussian quadrature nodes.

  Note: This function does NOT cacluate the derivative of f.
  Must pass in the derivative manually.
*/
Vector FE_LoadVector2D(const FEM2D& fem, real2DFunction f, const int n_gq, const int xDerivativeOrder, const int yDerivativeOrder);

/*
  \returns the FE mass matrix for a function "a" using one FEM2D.

  \param a: Function in the inner products of the mass matrix.
  \param n_gq: Number of Gaussian quadrature nodes.

  Note: Full matrix representation is inefficient here.
  Better to use sparse matrix representation.
*/
Matrix FE_MassMatrix2D(const FEM2D& fem, real2DFunction a, const int n_gq, const int xDerivativeOrder, const int yDerivativeOrder);

/*
  \returns the FE mass matrix for a function "a" using one FEM2D.

  \param a: Function in the inner products of the mass matrix.
  \param n_gq: Number of Gaussian quadrature nodes.

  Note: Full matrix representation is inefficient here.
  Better to use sparse matrix representation.
*/
Matrix FE_MassMatrix2D(const FEM2D& fem,
                       real2DFunction a,
                       const int n_gq, 
                       const int xDerivativeOrder1, const int yDerivativeOrder1, 
                       const int xDerivativeOrder2, const int yDerivativeOrder2);

/*
  Performs and L2 projection on FEM2D for a function f.

  \param n_gq: Number of Gaussian quadrature nodes.
*/
void L2_Projection2D(FEM2D& fem, real2DFunction f, const int n_gq);