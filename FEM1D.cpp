#include "Precompilied.h"
#include "FEM1D.h"
#include "LagrangeShapeFunctions1D.h"

FEM1D::FEM1D(const Mesh1D& FEmesh, int order)
  : mesh(FEmesh),
  meshSize(mesh.size),
  polynomialOrder(order),
  Ng(mesh.size* order + 1),
  Nu(Ng - mesh.numBoundaryNodes)
{
  const int& p = polynomialOrder;

  // Debug
  ASSERT(p > 0, "Polynomial order must be positive");
  ASSERT(mesh.numBoundaryNodes >= 0, "Boundary conditions have not been set up in mesh");

  // Create arrays
  FENodes = new FENode1D[Ng];
  connectivityMatrix.resize(meshSize);
  for (int i = 0; i < meshSize; ++i)
    connectivityMatrix[i].resize(p + 1);

  // Initialize arrays
  for (int elem = 0; elem < meshSize; ++elem)
  {
    // Initialize node coordinates
    const real xL = mesh(elem, EdgeType::Left).x;
    const real xR = mesh(elem, EdgeType::Right).x;
    for (int i = 0; i < p; ++i)
    {
      FENodes[elem * p + i].x = xL + i * (xR - xL) / p;
      connectivityMatrix[elem][i] = elem * p + i;
    }
    connectivityMatrix[elem][p] = (elem + 1) * p;

    // Initialize boundary conditions
    FENodes[elem * p].BC = mesh(elem, EdgeType::Left).BC;
    if ((int)mesh(elem, EdgeType::Left).BC >= 0)
      boundaryIndices.emplace_back(elem * p);
  }
#pragma warning(suppress: 6386)
  FENodes[meshSize * p].x = mesh(meshSize - 1, EdgeType::Right).x;
  FENodes[meshSize * p].BC = mesh(meshSize - 1, EdgeType::Right).BC;
  if ((int)mesh(meshSize - 1, EdgeType::Right).BC >= 0)
    boundaryIndices.emplace_back(meshSize * p);
}

FEM1D::FEM1D(const Mesh1D& FEmesh, int order, realFunction initialCondition)
  : FEM1D(FEmesh, order)
{
  const auto& f = initialCondition;

  // Interpolate initialCondition
  for (int i = 0; i < Ng; ++i)
    FENodes[i].u = f(FENodes[i].x);
}

FEM1D::FEM1D(FEM1D&& other) noexcept
  : mesh(other.mesh),
  meshSize(other.meshSize),
  polynomialOrder(other.polynomialOrder),
  Ng(other.Ng),
  Nu(other.Nu),
  boundaryIndices(std::move(other.boundaryIndices))
{
  FENodes = other.FENodes;
  connectivityMatrix = std::move(other.connectivityMatrix);

  other.FENodes = nullptr;
}

std::vector<int>& FEM1D::operator[](const int elementIndex)
{
  // Debug
  ASSERT(elementIndex >= 0, "Element index must be non-negative");
  ASSERT(elementIndex < meshSize, "Element index must be less than the number of elements");

  return connectivityMatrix[elementIndex];
}

const std::vector<int>& FEM1D::operator[](const int elementIndex) const
{
  // Debug
  ASSERT(elementIndex >= 0, "Element index must be non-negative");
  ASSERT(elementIndex < meshSize, "Element index must be less than the number of elements");

  return connectivityMatrix[elementIndex];
}

FENode1D& FEM1D::operator()(const int elementIndex, const int nodeIndex) const
{
  // Debug
  ASSERT(elementIndex >= 0, "Element index must be non-negative");
  ASSERT(elementIndex < meshSize, "Element index must be less than the number of elements");
  ASSERT(nodeIndex >= 0, "Node index must be non-negative");
  ASSERT(nodeIndex <= polynomialOrder, "Node index must be less than or equal to the polynomial order");

  return FENodes[connectivityMatrix[elementIndex][nodeIndex]];
}

FEM1D::~FEM1D()
{
  delete[] FENodes;
}

real FEM1D::evaluate(const real x) const
{
  return evaluate(x, 0);
}

real FEM1D::integrateElement(const int elementIndex, const int n_gq) const
{
  // Debug
  ASSERT(elementIndex >= 0, "Element index must be non-negative");
  ASSERT(elementIndex < meshSize, "Element index must be less than the number of elements");

  const int& K = elementIndex;
  std::vector<real> GLNodes = gauss1DNodesLocal(mesh, K, n_gq);
  std::vector<real> GLWeights = gauss1DWeightsLocal(mesh, K, n_gq);

  real sum = 0.0;
  for (int i = 0; i < n_gq; ++i)
  {
    real integrand = evaluate(GLNodes[i], elementIndex, 0);
    sum += GLWeights[i] * integrand;
  }
  return sum;
}

real FEM1D::integrate(const int n_gq) const
{
  real sum = 0;
  for (int K = 0; K < meshSize; ++K)
    sum += integrateElement(K, n_gq);
  return sum;
}

real FEM1D::evaluate(const real x, const int derivativeOrder) const
{
  // Search for element that x is an element of
  int K = -1;
  for (int i = 0; i < meshSize; ++i)
  {
    const real& xL = mesh(i, EdgeType::Left).x;
    const real& xR = mesh(i, EdgeType::Right).x;
    if (x >= xL && x <= xR)
    {
      K = i;
      break;
    }
  }
  ASSERT(K >= 0, "x must be in the domain of the mesh");

  // Evaluate at x
  return evaluate(x, K, derivativeOrder);
}

real FEM1D::evaluate(const real x, const int elementIndex, const int derivativeOrder) const
{
  const int& K = elementIndex;
  const real& xL = mesh(K, EdgeType::Left).x;
  const real& xR = mesh(K, EdgeType::Right).x;

  // Debug
  ASSERT(elementIndex >= 0, "Element index must be non-negative");
  ASSERT(elementIndex < meshSize, "Element index must be less than the number of elements");
  ASSERT(x >= xL && x <= xR, "x must be in the element at the specified elementIndex");

  real sum = 0.0;
  for (int j = 0; j < polynomialOrder + 1; ++j)
  {
    const int& i = connectivityMatrix[K][j];
    sum += FENodes[i].u * lagrangeShapeFunction1D(x, *this, K, j, derivativeOrder);
  }
  return sum;
}

void FEM1D::plot(const int n, const int derivativeOrder) const
{
  std::ofstream file("Out.txt");
  for (int K = 0; K < meshSize; ++K)
  {
    for (int i = 0; i < n; ++i)
    {
      const real& xL = mesh(K, EdgeType::Left).x;
      const real& xR = mesh(K, EdgeType::Right).x;
      const real x = xL + i * (xR - xL) / n;

      file << x << ", " << evaluate(x, K, derivativeOrder) << std::endl;
    }
  }
  file.close();
  std::cin.get();
  remove("Out.txt");
}
