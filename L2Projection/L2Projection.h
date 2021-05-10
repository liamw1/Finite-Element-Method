#pragma once
#include "Precompilied.h"
#include "LinearAlgebra/Matrix.h"
#include "Functions/Gauss-LegendreNodes.h"
#include "Functions/LagrangeShapeFunctions1D.h"
#include "Functions/LagrangeShapeFunctions2D.h"

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
  \returns the FE load vector for a function f.

  \param f: Function in the inner products of the load vector.
  \param n_gq: Number of Gaussian quadrature nodes.

  Note: This function does NOT cacluate the derivative of f.
  Must pass in the derivative manually.
*/
template<int N>
Vector* FE_LoadVector2D(const FEM2D<N>& fem, real2DFunction f, const int n_gq, const int xDerivativeOrder, const int yDerivativeOrder)
{
  const int& p = fem.polynomialOrder;

  Vector* b = new Vector(fem.Ng);
  for (int K = 0; K < fem.mesh.size; ++K)
  {
    std::vector<std::array<real, 2>> GLnodes = gauss2DNodesLocal(fem.mesh, K, n_gq);
    std::vector<real> GLweights = gauss2DWeightsLocal(fem.mesh, K, n_gq);

    for (int j = 0; j < (p + 1) * (p + 2) / 2; ++j)
    {
      // Calculate inner product between f and j-th shape function on K
      real innerProduct = 0.0;
      for (int i = 0; i < n_gq; ++i)
      {
        const real& x = GLnodes[i][0];
        const real& y = GLnodes[i][1];
        const real integrand = f(x, y) * lagrangeShapeFunction2D(x, y, fem, K, j, xDerivativeOrder, yDerivativeOrder);
        innerProduct += GLweights[i] * integrand;
      }
      (*b)[fem[K][j]] += innerProduct; // Accumulate to b
    }
  }
  return b;
}
template<int N>
Vector* FE_LoadVector2D(const FEM2D<N>& fem, std::function<real(real, real)> f, const int n_gq, const int xDerivativeOrder, const int yDerivativeOrder)
{
  const int& p = fem.polynomialOrder;

  Vector* b = new Vector(fem.Ng);
  for (int K = 0; K < fem.mesh.size; ++K)
  {
    std::vector<std::array<real, 2>> GLnodes = gauss2DNodesLocal(fem.mesh, K, n_gq);
    std::vector<real> GLweights = gauss2DWeightsLocal(fem.mesh, K, n_gq);

    for (int j = 0; j < (p + 1) * (p + 2) / 2; ++j)
    {
      // Calculate inner product between f and j-th shape function on K
      real innerProduct = 0.0;
      for (int i = 0; i < n_gq; ++i)
      {
        const real& x = GLnodes[i][0];
        const real& y = GLnodes[i][1];
        const real integrand = f(x, y) * lagrangeShapeFunction2D(x, y, fem, K, j, xDerivativeOrder, yDerivativeOrder);
        innerProduct += GLweights[i] * integrand;
      }
      (*b)[fem[K][j]] += innerProduct; // Accumulate to b
    }
  }
  return b;
}

/*
  \returns the FE mass matrix for a function "a" using one FEM2D.

  \param a: Function in the inner products of the mass matrix.
  \param n_gq: Number of Gaussian quadrature nodes.

  Note: Full matrix representation is inefficient here.
  Better to use sparse matrix representation.
*/
template<int N>
Matrix* FE_MassMatrix2D(const FEM2D<N>& fem, real2DFunction a, const int n_gq, const int xDerivativeOrder, const int yDerivativeOrder)
{
  return FE_MassMatrix2D(fem, fem, a, n_gq, xDerivativeOrder, yDerivativeOrder, xDerivativeOrder, yDerivativeOrder);
}
template<int N>
Matrix* FE_MassMatrix2D(const FEM2D<N>& fem,
  real2DFunction a,
  const int n_gq,
  const int xDerivativeOrder1, const int yDerivativeOrder1,
  const int xDerivativeOrder2, const int yDerivativeOrder2)
{
  return FE_MassMatrix2D(fem, fem, a, n_gq, xDerivativeOrder1, yDerivativeOrder1, xDerivativeOrder2, yDerivativeOrder2);
}
template<int N>
Matrix* FE_MassMatrix2D(const FEM2D<N>& fem, std::function<real(real, real)> a, const int n_gq, const int xDerivativeOrder, const int yDerivativeOrder)
{
  return FE_MassMatrix2D(fem, fem, a, n_gq, xDerivativeOrder, yDerivativeOrder, xDerivativeOrder, yDerivativeOrder);
}
template<int N>
Matrix* FE_MassMatrix2D(const FEM2D<N>& fem,
  std::function<real(real, real)> a,
  const int n_gq,
  const int xDerivativeOrder1, const int yDerivativeOrder1,
  const int xDerivativeOrder2, const int yDerivativeOrder2)
{
  return FE_MassMatrix2D(fem, fem, a, n_gq, xDerivativeOrder1, yDerivativeOrder1, xDerivativeOrder2, yDerivativeOrder2);
}

/*
  \returns the FE mass matrix for a function "a" using two FEM2Ds.

  \param a: Function in the inner products of the mass matrix.
  \param n_gq: Number of Gaussian quadrature nodes.

  Note: Full matrix representation is inefficient here.
  Better to use sparse matrix representation.
*/
template<int N, int M>
Matrix* FE_MassMatrix2D(const FEM2D<N>& fem1, const FEM2D<M>& fem2,
  real2DFunction a,
  const int n_gq,
  const int xDerivativeOrder, const int yDerivativeOrder)
{
  return FE_MassMatrix2D(fem1, fem2, a, n_gq, xDerivativeOrder, yDerivativeOrder, xDerivativeOrder, yDerivativeOrder);
}
template<int N, int M>
Matrix* FE_MassMatrix2D(const FEM2D<N>& fem1, const FEM2D<M>& fem2,
  real2DFunction a,
  const int n_gq,
  const int xDerivativeOrder1, const int yDerivativeOrder1,
  const int xDerivativeOrder2, const int yDerivativeOrder2)
{
  // Debug
  ASSERT(fem1.mesh.size == fem2.mesh.size, "FEM structures do not share the same mesh");
  ASSERT(fem1.mesh.numNodes == fem2.mesh.numNodes, "FEM structures do not share the same mesh");
  ASSERT(fem1.mesh.numEdges == fem2.mesh.numEdges, "FEM structures do not share the same mesh");

  const int& p1 = fem1.polynomialOrder;
  const int& p2 = fem2.polynomialOrder;

  Matrix* A = new Matrix(fem1.Ng, fem2.Ng);
  for (int K = 0; K < fem1.mesh.size; ++K)
  {
    std::vector<std::array<real, 2>> GLnodes = gauss2DNodesLocal(fem1.mesh, K, n_gq);
    std::vector<real> GLweights = gauss2DWeightsLocal(fem1.mesh, K, n_gq);

    for (int i = 0; i < (p1 + 1) * (p1 + 2) / 2; ++i)
      for (int j = 0; j < (p2 + 1) * (p2 + 2) / 2; ++j)
      {
        // Calculate the inner product between the i-th and j-th shape function on K
        real innerProduct = 0.0;
        for (int k = 0; k < n_gq; ++k)
        {
          const real& x = GLnodes[k][0];
          const real& y = GLnodes[k][1];
          const real integrand = a(x, y) * lagrangeShapeFunction2D(x, y, fem1, K, i, xDerivativeOrder1, yDerivativeOrder1)
            * lagrangeShapeFunction2D(x, y, fem2, K, j, xDerivativeOrder2, yDerivativeOrder2);
          innerProduct += GLweights[k] * integrand;
        }
        (*A)[fem1[K][i]][fem2[K][j]] += innerProduct; // Accumulate to M
      }
  }
  return A;
}
template<int N, int M>
Matrix* FE_MassMatrix2D(const FEM2D<N>& fem1, const FEM2D<M>& fem2,
  std::function<real(real, real)> a,
  const int n_gq,
  const int xDerivativeOrder, const int yDerivativeOrder)
{
  return FE_MassMatrix2D(fem1, fem2, a, n_gq, xDerivativeOrder, yDerivativeOrder, xDerivativeOrder, yDerivativeOrder);
}
template<int N, int M>
Matrix* FE_MassMatrix2D(const FEM2D<N>& fem1, const FEM2D<M>& fem2,
  std::function<real(real, real)> a,
  const int n_gq,
  const int xDerivativeOrder1, const int yDerivativeOrder1,
  const int xDerivativeOrder2, const int yDerivativeOrder2)
{
  // Debug
  ASSERT(fem1.mesh.size == fem2.mesh.size, "FEM structures do not share the same mesh");
  ASSERT(fem1.mesh.numNodes == fem2.mesh.numNodes, "FEM structures do not share the same mesh");
  ASSERT(fem1.mesh.numEdges == fem2.mesh.numEdges, "FEM structures do not share the same mesh");

  const int& p1 = fem1.polynomialOrder;
  const int& p2 = fem2.polynomialOrder;

  Matrix* A = new Matrix(fem1.Ng, fem2.Ng);
  for (int K = 0; K < fem1.mesh.size; ++K)
  {
    std::vector<std::array<real, 2>> GLnodes = gauss2DNodesLocal(fem1.mesh, K, n_gq);
    std::vector<real> GLweights = gauss2DWeightsLocal(fem1.mesh, K, n_gq);

    for (int i = 0; i < (p1 + 1) * (p1 + 2) / 2; ++i)
      for (int j = 0; j < (p2 + 1) * (p2 + 2) / 2; ++j)
      {
        // Calculate the inner product between the i-th and j-th shape function on K
        real innerProduct = 0.0;
        for (int k = 0; k < n_gq; ++k)
        {
          const real& x = GLnodes[k][0];
          const real& y = GLnodes[k][1];
          const real integrand = a(x, y) * lagrangeShapeFunction2D(x, y, fem1, K, i, xDerivativeOrder1, yDerivativeOrder1)
            * lagrangeShapeFunction2D(x, y, fem2, K, j, xDerivativeOrder2, yDerivativeOrder2);
          innerProduct += GLweights[k] * integrand;
        }
        (*A)[fem1[K][i]][fem2[K][j]] += innerProduct; // Accumulate to M
      }
  }
  return A;
}

/*
  Performs and L2 projection on FEM2D for a function f.

  \param n_gq: Number of Gaussian quadrature nodes.
*/
template<int N>
void L2_Projection2D(FEM2D<N>& fem, real2DFunction f, const int n_gq)
{
  const int u = 0;
  const int& p = fem.polynomialOrder;
  real2DFunction identityFunction = [](real x, real y) { return 1.0; };

  Matrix* M = FE_MassMatrix2D(fem, identityFunction, n_gq, 0, 0);
  Vector* b = FE_LoadVector2D(fem, f, n_gq, 0, 0);
  Vector coefficients = solve(*M, *b);

  for (int K = 0; K < fem.mesh.size; ++K)
    for (int j = 0; j < (p + 1) * (p + 2) / 2; ++j)
      fem(K, j)[u] = coefficients[fem[K][j]];
}