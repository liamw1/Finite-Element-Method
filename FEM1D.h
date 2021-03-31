#pragma once
#include "Precompilied.h"
#include "Mesh1D.h"
#include "Vector.h"
#include "Gauss-LegendreNodes.h"

class FEM1D
{
public:
  const Mesh1D& mesh;
  const int meshSize;
  const int polynomialOrder;
  const int Ng; // Number of FE nodes
  const int Nu; // Number of non-boundary FE nodes
  std::vector<int> boundaryIndices{};  // Indices of all boundary nodes, must be ordered!
  FENode1D* FENodes; // Stores all finite element nodes

  FEM1D() = delete;

  /*
    Generates a finite element method using a given mesh and a
    polynomial order.
  */
  FEM1D(const Mesh1D& FEmesh, const int order);

  /*
    Generates a finite element method using a given mesh, a
    polynomial order, and a (continuous) initial condition function.

    Will use Lagrange interpolation to approximate the initial condition function.
  */
  FEM1D(const Mesh1D& mesh, const int order, realFunction initialCondition);

  /*
    Move constructor.
  */
  FEM1D(FEM1D&& other) noexcept;

  FEM1D(const FEM1D& other) = delete;

  FEM1D& operator=(const FEM1D& other) = delete;

  /*
    \returns an array containing the indices of the FE nodes
    that belong to the element at the specified index.
  */
  std::vector<int>& operator[](const int elementIndex);

  /*
    \returns an array containing the indices of the FE nodes
    that belong to the element at the specified index.
  */
  const std::vector<int>& operator[](const int elementIndex) const;

  /*
    \returns the FE node at the specified indices.
    
    \param nodeIndex: Index of the nodes within a particular element.
  */
  FENode1D& operator()(const int elementIndex, const int nodeIndex) const;

  ~FEM1D();

  /*
    \returns the numerical approximation u_h(x).

    Note: This function is ineffiecient, as it searches for the element
    that x is in.  If the element is known, use other evaluate function.
  */
  real evaluate(const real x) const;

  /*
    \returns the numerical approximation u_h(x) for a given 
    derivative order.

    Note: This function is ineffiecient, as it searches for the element
    that x is in.  If the element is known, use other evaluate function.
  */
  real evaluate(const real x, const int derivativeOrder) const;
  
  /*
    \returns the numerical approximation u_h(x) for a given 
    derivative order.
    
    \param elementIndex: Index of element that contains x.
  */
  real evaluate(const real x, const int elementIndex, const int derivativeOrder) const;

  /*
    \returns the definite integral of the numerical approximation
    u_h(x) over the entire domain.
    
    \param n_gq: Number of Gaussian quadrature nodes.
  */
  real integrate(const int n_gq) const;

  /*
    Outputs a file called "Out.txt" that when graphed gives the
    d-th derivative of the FE interpolation.

    \param n: Number of points to plot per element.

    \param derivativeOrder: (optional)
  */
  void plot(const int n, const int derivativeOrder = 0) const;

private:
  /*
    A 2D array that stores which FE nodes belong to each mesh element.
    connectivityMatrix[i] gives an array containing
    the indices of the nodes that belong to the i-th element.
  */
  std::vector<std::vector<int>> connectivityMatrix;

  /*
    \returns the definite integral of the numerical approximation
    u_h(x) over a specified mesh element.
    
    \param n_gq: Number of Gaussian quadrature nodes.
  */
  real integrateElement(const int elementIndex, const int n_gq) const;
};