#pragma once
#include "Precompilied.h"
#include "Matrix.h"
#include "Gauss-LegendreNodes.h"
#include "LagrangeShapeFunctions1D.h"

/*
  \returns the FE load vector for an analytical function f.

  \param f: Analytical function to compare to.

  \param n_gq: Number of Gaussian quadrature nodes.

  \param derivativeOrder: Order of derivative on the Lagrange shape functions.  This function
  does NOT cacluate the derivative of f.  Must pass in the derivative manually.
*/
Vector FE_LoadVector1D(const FEM& fem, realFunction f, const int n_gq, const int derivativeOrder);

/*
  \returns the FE mass matrix using one FEM.

  \param n_gq: Number of Gaussian quadrature nodes.

  \param derivativeOrder: Order of derivative on the Lagrange shape functions.

  Note: Full matrix representation is inefficient here.  
  Better to use sparse matrix representation.
*/
Matrix FE_MassMatrix1D(const FEM& fem, const int n_gq, const int derivativeOrder);

/*
  \returns the FE mass matrix using one FEM.

  \param n_gq: Number of Gaussian quadrature nodes.

  \param derivativeOrder1: Order of derivative on the Lagrange shape functions on test functions v.

  \param derivativeOrder2: Order of derivative on the Lagrange shape functions on trial functions u.

  Note: Full matrix representation is inefficient here.
  Better to use sparse matrix representation.
*/
Matrix FE_MassMatrix1D(const FEM& fem, const int n_gq, const int derivativeOrder1, const int derivativeOrder2);

/*
  \returns the FE mass matrix for a function "a" using one FEM.

  \param a: Function in the inner products of the mass matrix.

  \param n_gq: Number of Gaussian quadrature nodes.

  \param derivativeOrder: Order of derivative on the Lagrange shape functions.

  Note: Full matrix representation is inefficient here.
  Better to use sparse matrix representation.
*/
Matrix FE_MassMatrix1D(const FEM& fem, realFunction a, const int n_gq, const int derivativeOrder);

/*
  \returns the FE mass matrix for a function "a" using one FEM.

  \param a: Function in the inner products of the mass matrix.

  \param n_gq: Number of Gaussian quadrature nodes.

  \param derivativeOrder1: Order of derivative on the Lagrange shape functions on test functions v.

  \param derivativeOrder2: Order of derivative on the Lagrange shape functions on trial functions u.

  Note: Full matrix representation is inefficient here.
  Better to use sparse matrix representation.
*/
Matrix FE_MassMatrix1D(const FEM& fem, realFunction a, const int n_gq, const int derivativeOrder1, const int derivativeOrder2);

/*
  Performs and L2 projection on fem for a function f.

  \param n_gq: Number of Gaussian quadrature nodes.

  \param derivativeOrder: Order of derivative on the Lagrange shape functions.
*/
void L2_Projection(FEM& fem, realFunction f, const int n_gq, const int derivativeOrder);