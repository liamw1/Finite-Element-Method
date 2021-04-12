#include "Precompilied.h"
#include "ErrorAnalysis.h"

real FE_AbsoluteError1D(const real x, const int elementIndex, const FEM1D& fem, real1DFunction f, const int derivativeOrder)
{
  const int& LEFT = 0, RIGHT = fem.polynomialOrder;
  const int& K = elementIndex;
  const real& xL = fem(K, LEFT).x;
  const real& xR = fem(K, RIGHT).x;
  ASSERT(x >= xL && x <= xR, "x must be in the element at the specified elementIndex");

  return abs(f(x) - fem.evaluate(x, elementIndex, derivativeOrder));
}

void plotAbsoluteError1D(const FEM1D& fem, real1DFunction f, const int n, const int derivativeOrder)
{
  const int& LEFT = 0, RIGHT = fem.polynomialOrder;

  std::ofstream file("Out.txt");
  for (int K = 0; K < fem.meshSize; ++K)
    for (int i = 0; i < n; ++i)
    {
      const real& xL = fem(K, LEFT).x;
      const real& xR = fem(K, RIGHT).x;
      const real x = xL + i * (xR - xL) / n;

      file << x << ", " << FE_AbsoluteError1D(x, K, fem, f, derivativeOrder) << std::endl;
    }
  file.close();
  std::cin.get();
  remove("Out.txt");
}

real FE_Error1DLocal(const int elementIndex, const FEM1D& fem, real1DFunction f, const int n_gq, const int derivativeOrder)
{
  const int& K = elementIndex;
  std::vector<real> GLnodes = gauss1DNodesLocal(fem.mesh, K, n_gq);
  std::vector<real> GLweights = gauss1DWeightsLocal(fem.mesh, K, n_gq);

  real sum = 0.0;
  for (int i = 0; i < n_gq; ++i)
  {
    const real err = FE_AbsoluteError1D(GLnodes[i], K, fem, f, derivativeOrder);
    sum += GLweights[i] * err * err;
  }

  return sqrtl(sum);
}

real FE_Error1DGlobal(const FEM1D& fem, real1DFunction f, const int n_gq, const int derivativeOrder)
{
  real sum = 0.0;
  for (int K = 0; K < fem.meshSize; ++K)
  {
    const real localNorm = FE_Error1DLocal(K, fem, f, n_gq, derivativeOrder);
    sum += localNorm * localNorm;
  }
  return sqrtl(sum);
}

real FE_AbsoluteError2D(const real x, const real y,
                        const int elementIndex,
                        const FEM2D& fem,
                        real2DFunction f,
                        const int xDerivativeOrder, const int yDerivativeOrder)
{
  const int& K = elementIndex;
  ASSERT(fem.isInTriangle(x, y, K), "x must be in the element at the specified elementIndex");

  return abs(f(x, y) - fem.evaluate(x, y, elementIndex, xDerivativeOrder, yDerivativeOrder));
}

void plotAbsoluteError2D(const FEM2D& fem, real2DFunction f, const int n, const int xDerivativeOrder, const int yDerivativeOrder)
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

real FE_Error2DLocal(const int elementIndex, const FEM2D& fem, real2DFunction f, const int n_gq, const int xDerivativeOrder, const int yDerivativeOrder)
{
  const int& K = elementIndex;
  std::vector<std::array<real, 2>> GLnodes = gauss2DNodesLocal(fem.mesh, K, n_gq);
  std::vector<real> GLweights = gauss2DWeightsLocal(fem.mesh, K, n_gq);

  real sum = 0.0;
  for (int i = 0; i < n_gq; ++i)
  {
    const real& x = GLnodes[i][0];
    const real& y = GLnodes[i][1];
    const real err = FE_AbsoluteError2D(x, y, K, fem, f, xDerivativeOrder, yDerivativeOrder);
    sum += GLweights[i] * err * err;
  }

  return sqrtl(sum);
}

real FE_Error2DGlobal(const FEM2D& fem, real2DFunction f, const int n_gq, const int xDerivativeOrder, const int yDerivativeOrder)
{
  real sum = 0.0;
  for (int K = 0; K < fem.mesh.size; ++K)
  {
    const real localNorm = FE_Error2DLocal(K, fem, f, n_gq, xDerivativeOrder, yDerivativeOrder);
    sum += localNorm * localNorm;
  }
  return sqrtl(sum);
}
