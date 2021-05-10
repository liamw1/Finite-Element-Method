#pragma once
#include "Precompilied.h"
#include "Meshing/1D/FEM1D.h"
#include "Meshing/2D/FEM2D.h"
#include "Functions/Gauss-LegendreNodes.h"

/*
  \returns absolute error between FE interpolation and analytical function f.

  \param elementIndex: Index of element that contains x.

  \param f: Analytical function to compare to.

  \param derivativeOrder: Order of derivative of FE interpolation.  This function
  does NOT cacluate the derivative of f.  Must pass in the derivative manually.
*/
real FE_AbsoluteError1D(const real x, const int elementIndex, const FEM1D& fem, real1DFunction f, const int derivativeOrder);

/*
  Outputs a file called "Out.txt" that when graphed gives the
  absolute error between FE interpolation and analytical function f.

  \param f: Analytical function to compare to.

  \param n: Number of points to plot per element.

  \param derivativeOrder: Order of derivative of FE interpolation.  This function
  does NOT cacluate the derivative of f.  Must pass in the derivative manually.
*/
void plotAbsoluteError1D(const FEM1D& fem, real1DFunction f, const int n, const int derivativeOrder);

/*
  \returns the L2 norm of the absolute error over a given element.  Integration done
  using Gaussian quadrature.

  \param f: Analytical function to compare to.

  \param n_gq: Number of Gaussian quadrature nodes.

  \param derivativeOrder: Order of derivative of FE interpolation.  This function
  does NOT cacluate the derivative of f.  Must pass in the derivative manually.
*/
real FE_Error1DLocal(const int elementIndex, const FEM1D& fem, real1DFunction f, const int n_gq, const int derivativeOrder);

/*
  \returns the L2 norm of the absolute error over the entire domain.
  Integration done using Gaussian quadrature.

  \param f: Analytical function to compare to.

  \param n_gq: Number of Gaussian quadrature nodes.

  \param derivativeOrder: Order of derivative of FE interpolation.  This function
  does NOT cacluate the derivative of f.  Must pass in the derivative manually.
*/
real FE_Error1DGlobal(const FEM1D& fem, real1DFunction f, const int n_gq, const int derivativeOrder);



/*
  \returns absolute error between FE interpolation and analytical function f.

  \param elementIndex: Index of element that contains x.

  \param f: Analytical function to compare to.

  Note: This function does NOT cacluate the derivatives of f.  Must pass in the derivative manually.
*/
template<int N>
real FE_AbsoluteError2D(const int varIndex,
                        const real x, const real y,
                        const int elementIndex,
                        const FEM2D<N>& fem,
                        real2DFunction f,
                        const int xDerivativeOrder, const int yDerivativeOrder)
{
  const int& K = elementIndex;
  ASSERT(fem.isInTriangle(x, y, K), "x must be in the element at the specified elementIndex");

  return abs(f(x, y) - fem.evaluate(varIndex, x, y, elementIndex, xDerivativeOrder, yDerivativeOrder));
}

/*
  Outputs a file called "Plot.txt" that when graphed gives the
  absolute error between FE interpolation and analytical function f.

  \param f: Analytical function to compare to.

  \param n: Number of points to plot per direction per element.

  Note: This function does NOT cacluate the derivatives of f.  Must pass in the derivative manually.
*/
template<int N>
void plotAbsoluteError2D(const FEM2D<N>& fem, real2DFunction f, const int n, const int xDerivativeOrder, const int yDerivativeOrder)
{
  std::ofstream file("Plot.txt");

  for (int K = 0; K < fem.mesh.size; ++K)
  {
    const MeshNode2D& A1 = fem.mesh(K, 0);
    const MeshNode2D& A2 = fem.mesh(K, 1);
    const MeshNode2D& A3 = fem.mesh(K, 2);

    // Create transformation matrix
    Matrix B = Matrix(2);
    B[0][0] = A2.x - A1.x;
    B[0][1] = A3.x - A1.x;
    B[1][0] = A2.y - A1.y;
    B[1][1] = A3.y - A1.y;

    for (int i = 0; i <= n; ++i)
      for (int j = 0; j <= n - i; ++j)
      {
        // Create reference points
        const real tx = (real)i / n;
        const real ty = (real)j / n;

        // Transform points from reference domain to element K
        const real x = B[0][0] * tx + B[0][1] * ty + A1.x;
        const real y = B[1][0] * tx + B[1][1] * ty + A1.y;

        file << x << ", " << y << ", " << FE_AbsoluteError2D(x, y, K, fem, f, xDerivativeOrder, yDerivativeOrder) << std::endl;
      }
  }
  file.close();
  std::cout << "Data written to \"Plot.txt\".  Press any key to continue" << std::endl;
  std::cin.get();
  remove("Plot.txt");
}

/*
  \returns the L2 norm of the absolute error over a given element.  Integration done
  using Gaussian quadrature.

  \param f: Analytical function to compare to.

  \param n_gq: Number of Gaussian quadrature nodes.

  Note: This function does NOT cacluate the derivatives of f.  Must pass in the derivative manually.
*/
template<int N>
real FE_Error2DLocal(const int varIndex, const int elementIndex, const FEM2D<N>& fem, real2DFunction f, const int n_gq, const int xDerivativeOrder, const int yDerivativeOrder)
{
  const int& K = elementIndex;
  std::vector<std::array<real, 2>> GLnodes = gauss2DNodesLocal(fem.mesh, K, n_gq);
  std::vector<real> GLweights = gauss2DWeightsLocal(fem.mesh, K, n_gq);

  real sum = 0.0;
  for (int i = 0; i < n_gq; ++i)
  {
    const real& x = GLnodes[i][0];
    const real& y = GLnodes[i][1];
    const real err = FE_AbsoluteError2D(varIndex, x, y, K, fem, f, xDerivativeOrder, yDerivativeOrder);
    sum += GLweights[i] * err * err;
  }

  return sqrtl(sum);
}

/*
  \returns the L2 norm of the absolute error over the entire domain.
  Integration done using Gaussian quadrature.

  \param f: Analytical function to compare to.

  \param n_gq: Number of Gaussian quadrature nodes.

  Note: This function does NOT cacluate the derivatives of f.  Must pass in the derivative manually.
*/
template<int N>
real FE_Error2DGlobal(const int varIndex, const FEM2D<N>& fem, real2DFunction f, const int n_gq, const int xDerivativeOrder, const int yDerivativeOrder)
{
  real sum = 0.0;
  for (int K = 0; K < fem.mesh.size; ++K)
  {
    const real localNorm = FE_Error2DLocal(varIndex, K, fem, f, n_gq, xDerivativeOrder, yDerivativeOrder);
    sum += localNorm * localNorm;
  }
  return sqrtl(sum);
}