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