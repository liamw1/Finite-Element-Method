#include "Precompilied.h"
#include "FEM2D.h"
#include "LagrangeShapeFunctions2D.h"

FEM2D::FEM2D(const Mesh2D& FEmesh, const int order)
  : mesh(FEmesh),
    polynomialOrder(order),
    Ng(mesh.numNodes + (order - 1) * mesh.numEdges + (order - 1) * (order - 2) * mesh.size / 2),
    connectivityMatrix(mesh.size, (order + 2) * (order + 1) / 2)
{
  const int& p = polynomialOrder;

  // Debug
  ASSERT(p > 0, "Polynomial order must be positive");

  FENodes = new FENode2D[Ng];

  // Handle vertices first
  for (int n = 0; n < mesh.numNodes; ++n)
  {
#pragma warning(suppress: 6386)
    FENodes[n].x = mesh.meshNodes[n].x;
    FENodes[n].y = mesh.meshNodes[n].y;
    FENodes[n].BC = mesh.meshNodes[n].BC;
  }
  for (int K = 0; K < mesh.size; ++K)
    for (int v = 0; v < 3; ++v)
      connectivityMatrix[K][v] = mesh.connectivityMatrix[K][v];

  // Then handle nodes along edges
  for (int e = 0; e < mesh.numEdges; ++e)
  {
    // Grab nodes that form edge
    const MeshNode2D& A1 = mesh.meshNodes[mesh.edgeArray[e][0]];
    const MeshNode2D& A2 = mesh.meshNodes[mesh.edgeArray[e][1]];

    // Create p-1 equally spaced nodes along the line passing through A1,A2
    for (int i = 0; i < p - 1; ++i)
    {
      const real t = (i + 1.0) / p;
#pragma warning(suppress: 6386)
      FENodes[mesh.numNodes + e * (p - 1) + i].x = (1 - t) * A1.x + t * A2.x;
      FENodes[mesh.numNodes + e * (p - 1) + i].y = (1 - t) * A1.y + t * A2.y;
    }
  }
  for (int K = 0; K < mesh.size; ++K)
  {
    int e = 0;
    for (int i = 1; i < 3; ++i)
      for (int j = i - 1; j >= 0; --j)
      {
        // Grab indices of potential edge-forming nodes
        const int& I = mesh.connectivityMatrix[K][i];
        const int& J = mesh.connectivityMatrix[K][j];

        if (mesh.edgeTypeMatrix[I][J] > 0) // Check if they form edge
        {
          const int E = mesh.edgeMatrix[I][J] - 1;
          ASSERT(E + 1 > 0, "Nodes I and J do not form an edge");

          for (int l = 0; l < p - 1; ++l)
            connectivityMatrix[K][3 + e * (p - 1) + l] = mesh.numNodes + E * (p - 1) + l;
          ++e;
        }
      }
  }

  // Finally, handle nodes on the interior
  int nodeIndex = mesh.numNodes + mesh.numEdges * (p - 1);
  for (int K = 0; K < mesh.size; ++K)
  {
    const MeshNode2D& A1 = mesh(K, 0);
    const MeshNode2D& A2 = mesh(K, 1);
    const MeshNode2D& A3 = mesh(K, 2);

    // Create transformation matrix from reference domain to K
    Matrix B = Matrix(2);
    B[0][0] = A2.x - A1.x;
    B[0][1] = A3.x - A1.x;
    B[1][0] = A2.y - A1.y;
    B[1][1] = A3.y - A1.y;

    int localNodeIndex = 3 * p;
    for (int i = 0; i < p - 2; ++i)
      for (int j = 0; j < p - 2 - i; ++j)
      {
        // Create reference points
        const real tx = (real)(i + 1) / p;
        const real ty = (real)(j + 1) / p;

        // Transform points from reference domain to element K
        FENodes[nodeIndex].x = B[0][0] * tx + B[0][1] * ty + A1.x;
        FENodes[nodeIndex].y = B[1][0] * tx + B[1][1] * ty + A1.y;

        connectivityMatrix[K][localNodeIndex] = nodeIndex;
        ++nodeIndex;
        ++localNodeIndex;
      }
  }
}

FEM2D::FEM2D(const Mesh2D& FEmesh, const int order, real2DFunction initialCondition)
  : FEM2D(FEmesh, order)
{
  const auto& f = initialCondition;

  // Interpolate initialCondition
  for (int i = 0; i < Ng; ++i)
    FENodes[i].u = f(FENodes[i].x, FENodes[i].y);
}

FEM2D::FEM2D(FEM2D&& other) noexcept
  : mesh(other.mesh),
    polynomialOrder(other.polynomialOrder),
    Ng(other.Ng),
    connectivityMatrix(std::move(connectivityMatrix))
{
  FENodes = other.FENodes;
  other.FENodes = nullptr;
}

const Container<int>& FEM2D::operator[](const int elementIndex) const
{
  // Debug
  ASSERT(elementIndex >= 0, "Element index must be non-negative");
  ASSERT(elementIndex < mesh.size, "Element index must be less than the number of elements");

  return connectivityMatrix[elementIndex];
}

FENode2D& FEM2D::operator()(const int elementIndex, const int nodeIndex) const
{
  // Debug
  ASSERT(elementIndex >= 0, "Element index must be non-negative");
  ASSERT(elementIndex < mesh.size, "Element index must be less than the number of elements");
  ASSERT(nodeIndex >= 0, "Node index must be non-negative");
  ASSERT(nodeIndex < (polynomialOrder + 1) * (polynomialOrder + 2) / 2, "Node index must be less than the number of nodes per element");

  return FENodes[connectivityMatrix[elementIndex][nodeIndex]];
}

real FEM2D::evaluate(const real x, const real y) const
{
  return evaluate(x, y, 0, 0);
}

real FEM2D::evaluate(const real x, const real y, const int xDerivativeOrder, const int yDerivativeOrder) const
{
  // Search for element that x is an element of
  int K = -1;
  for (int i = 0; i < mesh.size; ++i)
    if (isInTriangle(x, y, i))
    {
      K = i;
      break;
    }
  ASSERT(K >= 0, "x must be in the domain of the mesh");

  // Evaluate at x
  return evaluate(x, y, K, xDerivativeOrder, yDerivativeOrder);
}

real FEM2D::evaluate(const real x, const real y, const int elementIndex, const int xDerivativeOrder, const int yDerivativeOrder) const
{
  const int& K = elementIndex;
  const int& p = polynomialOrder;

  // Debug
  ASSERT(K >= 0, "Element index must be non-negative");
  ASSERT(K < mesh.size, "Element index must be less than the number of elements");
  ASSERT(isInTriangle(x, y, K), "x must be in the element at the specified elementIndex");

  real sum = 0.0;
  for (int j = 0; j < (p + 1) * (p + 2) / 2; ++j)
  {
    const int& i = connectivityMatrix[K][j];
    sum += FENodes[i].u * lagrangeShapeFunction2D(x, y, *this, K, j, xDerivativeOrder, yDerivativeOrder);
  }
  return sum;
}

void FEM2D::plot(const int n) const
{
  plot(n, 0, 0);
}

void FEM2D::plot(const int n, const int xDerivativeOrder, const int yDerivativeOrder) const
{
  std::ofstream file("Plot.txt");

  for (int K = 0; K < mesh.size; ++K)
  {
    const MeshNode2D& A1 = mesh(K, 0);
    const MeshNode2D& A2 = mesh(K, 1);
    const MeshNode2D& A3 = mesh(K, 2);

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

        file << x << ", " << y << ", " << evaluate(x, y, K, xDerivativeOrder, yDerivativeOrder) << std::endl;
      }
  }
  file.close();
  std::cout << "Data written to \"Plot.txt\".  Press any key to continue" << std::endl;
  std::cin.get();
  remove("Plot.txt");
}

bool FEM2D::isInTriangle(const real x, const real y, const int elementIndex) const
{
  const int& K = elementIndex;
  const MeshNode2D& A1 = mesh(K, 0);
  const MeshNode2D& A2 = mesh(K, 1);
  const MeshNode2D& A3 = mesh(K, 2);

  // Create transformation matrix from reference domain to K
  Matrix B = Matrix(2);
  B[0][0] = A2.x - A1.x;  B[0][1] = A3.x - A1.x;
  B[1][0] = A2.y - A1.y;  B[1][1] = A3.y - A1.y;

  // Invert matrix
  const real determinant = B[0][0] * B[1][1] - B[0][1] * B[1][0];
  ASSERT(determinant != 0.0, "Transformation matrix is singular");
  real temp = B[0][0];
  B[0][0] = B[1][1];
  B[1][1] = temp;
  B[0][1] *= -1;
  B[1][0] *= -1;
  B /= determinant;

  const real tx = B[0][0] * (x - A1.x) + B[0][1] * (y - A1.y);
  const real ty = B[1][0] * (x - A1.x) + B[1][1] * (y - A1.y);

  if (tx > -1.0 * TOLERANCE && ty > -1.0 * TOLERANCE && ty < 1.0 - tx + TOLERANCE)
    return true;
  else
    return false;
}