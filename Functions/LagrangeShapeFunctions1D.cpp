#include "Precompilied.h"
#include "LagrangeShapeFunctions1D.h"

real lagrangeShapeFunction1D(const real x,
                             const FEM1D& fem,
                             const int elementIndex,
                             const int nodeIndex,
                             const int derivativeOrder)
{
  const int& p = fem.polynomialOrder;
  const int& K = elementIndex;
  const int& j = nodeIndex;
  const real& xL = fem(K, 0).x;
  const real& xR = fem(K, p).x;
  const real scaling = 2.0 / (xR - xL);
  const real t = scaling * (x - xL) - 1.0;

  ASSERT(x >= xL && x <= xR, "x must be in the element at the specified elementIndex");

  // Mapping each node coordinate onto [-1, 1]
  std::vector<real> refNodes = std::vector<real>(p + 1);
  for (int i = 0; i < p + 1; ++i)
    refNodes[i] = scaling * (fem(K, i).x - xL) - 1.0;

  // Store the j-th node and remove it from list so that it's skipped over later
  const real t_j = refNodes[j];
  refNodes.erase(refNodes.begin() + j);

  ASSERT(t_j >= -1 && t_j <= 1, "The reference node t_j must be in the interval [-1, 1]");

  // A factor of "scaling" is acquired with each derivative
  real multiplier = 1.0;
  for (int i = 0; i < derivativeOrder; ++i)
    multiplier *= scaling;

  return multiplier * refLagrangePolynomial1D(t, t_j, refNodes, derivativeOrder);
}

real refLagrangePolynomial1D(const real t, const real t_j, std::vector<real>& refNodes, const int derivativeOrder)
{
  if (derivativeOrder == 0)
  {
    // Calculate Lagrange product
    real prod = 1.0;
    for (unsigned int i = 0; i < refNodes.size(); ++i)
    {
      const real& t_i = refNodes[i];
      prod *= (t - t_i) / (t_j - t_i);
    }
    return prod;
  }
  else
  {
    real sum = 0.0;
    for (unsigned int i = 0; i < refNodes.size(); ++i)
    {
      const real& t_i = refNodes[i];

      // Create copy of list with i-th node removed so that it's skipped over in the inner sums
      std::vector<real> refNodesCopy = refNodes;
      refNodesCopy.erase(refNodesCopy.begin() + i);

      // Recursively calculate inner sums
      sum += refLagrangePolynomial1D(t, t_j, refNodesCopy, derivativeOrder - 1) / (t_j - t_i);
    }
    return sum;
  }
}