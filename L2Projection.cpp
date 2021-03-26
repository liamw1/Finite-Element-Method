#include "Precompilied.h"
#include "L2Projection.h"

static real identityFunction(real x) { return 1.0; }

Vector FE_LoadVector1D(const FEM1D& fem, realFunction f, const int n_gq, const int derivativeOrder)
{
  Vector b = Vector(fem.Ng);
  for (int K = 0; K < fem.meshSize; ++K)
  {
    std::vector<real> GLnodes = gauss1DNodesLocal(fem.mesh, K, n_gq);
    std::vector<real> GLweights = gauss1DWeightsLocal(fem.mesh, K, n_gq);

    for (int j = 0; j < fem.polynomialOrder + 1; ++j)
    {
      // Calculate inner product between f and j-th shape function on K
      real innerProduct = 0.0;
      for (int i = 0; i < n_gq; ++i)
      {
        const real x = GLnodes[i];
        const real integrand = f(x) * lagrangeShapeFunction1D(x, fem, K, j, derivativeOrder);
        innerProduct += GLweights[i] * integrand;
      }
      b[fem[K][j]] += innerProduct; // Accumulate to b
    }
  }
  return b;
}

Matrix FE_MassMatrix1D(const FEM1D& fem, const int n_gq, const int derivativeOrder)
{
  return FE_MassMatrix1D(fem, identityFunction, n_gq, derivativeOrder, derivativeOrder);
}

Matrix FE_MassMatrix1D(const FEM1D& fem, const int n_gq, const int derivativeOrder1, const int derivativeOrder2)
{
  return FE_MassMatrix1D(fem, identityFunction, n_gq, derivativeOrder1, derivativeOrder2);
}

Matrix FE_MassMatrix1D(const FEM1D& fem, realFunction a, const int n_gq, const int derivativeOrder)
{
  return FE_MassMatrix1D(fem, a, n_gq, derivativeOrder, derivativeOrder);
}

Matrix FE_MassMatrix1D(const FEM1D& fem, realFunction a, const int n_gq, const int derivativeOrder1, const int derivativeOrder2)
{
  Matrix M = Matrix(fem.Ng);
  for (int K = 0; K < fem.meshSize; ++K)
  {
    std::vector<real> GLnodes = gauss1DNodesLocal(fem.mesh, K, n_gq);
    std::vector<real> GLweights = gauss1DWeightsLocal(fem.mesh, K, n_gq);

    for (int i = 0; i < fem.polynomialOrder + 1; ++i)
      for (int j = 0; j < fem.polynomialOrder + 1; ++j)
      {
        // Calculate the inner product between the i-th and j-th shape function on K
        real innerProduct = 0.0;
        for (int k = 0; k < n_gq; ++k)
        {
          const real x = GLnodes[k];
          const real integrand = a(x) * lagrangeShapeFunction1D(x, fem, K, i, derivativeOrder1)
            * lagrangeShapeFunction1D(x, fem, K, j, derivativeOrder2);
          innerProduct += GLweights[k] * integrand;
        }
        M[fem[K][i]][fem[K][j]] += innerProduct; // Accumulate to M
      }
  }
  return M;
}

void L2_Projection(FEM1D& fem, realFunction f, const int n_gq, const int derivativeOrder)
{
  Matrix M = FE_MassMatrix1D(fem, n_gq, derivativeOrder);
  Vector b = FE_LoadVector1D(fem, f, n_gq, derivativeOrder);
  Vector coefficients = solve(M, b);

  for (int K = 0; K < fem.meshSize; ++K)
    for (int j = 0; j < fem.polynomialOrder + 1; ++j)
      fem(K, j).u = coefficients[fem[K][j]];
}
