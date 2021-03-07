#pragma once
#include "Precompilied.h"
#include "FEM.h"
#include "Gauss-LegendreNodes.h"

/*
  \returns absolute error between FE interpolation and analytical function f.

  \param elementIndex: Index of element that contains x.

  \param f: Analytical function to compare to.

  \param derivativeOrder: Order of derivative of FE interpolation.  This function
  does NOT cacluate the derivative of f.  Must pass in the derivative manually.
*/
real FE_AbsoluteError1D(const real x, const int elementIndex, const FEM& fem, realFunction f, const int derivativeOrder);

/*
  Outputs a file called "Out.txt" that when graphed gives the
  absolute error between FE interpolation and analytical function f.

  \param f: Analytical function to compare to.

  \param n: Number of points to plot per element.

  \param derivativeOrder: Order of derivative of FE interpolation.  This function
  does NOT cacluate the derivative of f.  Must pass in the derivative manually.
*/
void plotAbsoluteError1D(const FEM& fem, realFunction f, const int n, const int derivativeOrder);

/*
  \returns the L2 norm of the absolute error over a given element.  Integration done
  using Gaussian quadrature.

  \param f: Analytical function to compare to.

  \param n_gq: Number of Gaussian quadrature nodes.

  \param derivativeOrder: Order of derivative of FE interpolation.  This function
  does NOT cacluate the derivative of f.  Must pass in the derivative manually.
*/
real FE_Error1DLocal(const int elementIndex, const FEM& fem, realFunction f, const int n_gq, const int derivativeOrder);

/*
  \returns the L2 norm of the absolute error over the entire domain.
  Integration done using Gaussian quadrature.

  \param f: Analytical function to compare to.

  \param n_gq: Number of Gaussian quadrature nodes.

  \param derivativeOrder: Order of derivative of FE interpolation.  This function
  does NOT cacluate the derivative of f.  Must pass in the derivative manually.
*/
real FE_Error1DGlobal(const FEM& fem, realFunction f, const int n_gq, const int derivativeOrder);