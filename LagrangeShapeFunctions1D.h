#pragma once
#include "Precompilied.h"
#include "FEM1D.h"

/*
  Calculates a 1D Lagrange shape function for an element K 
  using equally spaced local nodes.

  \param element: The index of an element within fem.
  \param node: The index of a node within the given element.
*/
real lagrangeShapeFunction1D(const real x, const FEM1D& fem, const int element, const int node, const int derivativeOrder);

/*
  Recursively calculates the nth derivative of the j-th Lagrange
  basis polynomial for p + 1 equally spaced local nodes on [-1, 1].
  Formula found here: https://en.wikipedia.org/wiki/Lagrange_polynomial#Derivatives

  \param t: The result when x is mapped onto [-1, 1].
  \param t_j: The coordinate of the j-th node when mapped to [-1, 1].
  \param refNodes: An array of node coordinates that have been mapped to [-1, 1].
*/
real refLagrangePolynomial1D(const real t, const real t_j, std::vector<real>& refNodes, const int derivativeOrder);
