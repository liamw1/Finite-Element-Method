#pragma once
#include "Precompilied.h"
#include "Mesh2D.h"
#include "Matrix.h"

/*
  2D finite element structure.

  Assumes triangular mesh.
*/
class FEM2D
{
public:
  const Mesh2D& mesh;
  const int polynomialOrder;
  const int Ng;
  FENode2D* FENodes;

  FEM2D() = delete;

  /*
    Generates a finite element method using a given mesh and a
    polynomial order.
  */
  FEM2D(const Mesh2D& FEmesh, const int order);

  /*
    Generates a finite element method using a given mesh, a
    polynomial order, and a (continuous) initial condition function.

    Will use Lagrange interpolation to approximate the initial condition function.
  */
  FEM2D(const Mesh2D& FEmesh, const int order, real2DFunction initialCondition);

  /*
    Move constructor.
  */
  FEM2D(FEM2D&& other) noexcept;

  FEM2D(const FEM2D& other) = delete;

  FEM2D& operator=(const FEM2D& other) = delete;

  const Container<int>& operator[](const int elementIndex) const;

  /*
    \returns the numerical approximation u_h(x, y).

    Note: This function is ineffiecient, as it searches for the element
    that (x, y) belongs to.  If the element is known, use other evaluate function.
  */
  real evaluate(const real x, const real y) const;
  real evaluate(const real x, const real y, const int xDerivativeOrder, const int yDerivativeOrder) const;

  /*
    \returns the numerical approximation u_h(x, y) for a given
    derivative order.

    \param elementIndex: Index of element that contains (x, y).
  */
  real evaluate(const real x, const real y, const int elementIndex, const int xDerivativeOrder, const int yDerivativeOrder) const;

  /*
    Outputs a file called "Plot.txt" that when graphed gives the
    d-th derivative of the FE interpolation.

    \param n: Number of points to plot per direction per element.
  */
  void plot(int n) const;
  void plot(int n, const int xDerivativeOrder, const int yDerivativeOrder) const;

  /*
    Determines if the given point (x, y) is inside the specified element.
  */
  bool isInTriangle(const real x, const real y, const int elementIndex) const;

private:
  /*
    A 2D array that stores which FE nodes belong to each mesh element.
    connectivityMatrix[i] gives an array containing
    the indices of the nodes that belong to the i-th element.
  */
  Array2D<int> connectivityMatrix;
};