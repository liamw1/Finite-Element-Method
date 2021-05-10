#pragma once
#include "Precompilied.h"
#include "Mesh2D.h"
#include "LinearAlgebra/Matrix.h"

/*
  2D finite element structure.

  Assumes triangular mesh.
*/
template <int N>
class FEM2D
{
public:
  const Mesh2D& mesh;
  const int polynomialOrder;
  const int Ng;                        // Number of FE nodes
  std::vector<int> boundaryIndices{};  // Indices of all boundary nodes, must be ordered!
  FENode2D<N>* FENodes;

  FEM2D() = delete;

  /*
    Generates a finite element method using a given mesh and a
    polynomial order.
  */
  FEM2D(const Mesh2D& FEmesh, const int order)
    : mesh(FEmesh),
    polynomialOrder(order),
    Ng(mesh.numNodes + (order - 1) * mesh.numEdges + (order - 1) * (order - 2) * mesh.size / 2),
    connectivityMatrix(mesh.size, (order + 2)* (order + 1) / 2)
  {
    const int& p = polynomialOrder;

    // Debug
    ASSERT(N > 0, "Number of variables must be positive");
    ASSERT(p > 0, "Polynomial order must be positive");
    ASSERT(mesh.numBoundaryNodes >= 0, "Boundary conditions have not been set up in mesh");

    FENodes = new FENode2D<N>[Ng];

    // Handle vertices first
    for (int n = 0; n < mesh.numNodes; ++n)
    {
#pragma warning(suppress: 6386)
      FENodes[n].x = mesh.meshNodes[n].x;
      FENodes[n].y = mesh.meshNodes[n].y;
      FENodes[n].BC = mesh.meshNodes[n].BC;
      FENodes[n].isCorner = mesh.meshNodes[n].isCorner;
      if ((int)FENodes[n].BC >= 0)
        boundaryIndices.emplace_back(n);
    }
    for (int K = 0; K < mesh.size; ++K)
      for (int v = 0; v < 3; ++v)
        connectivityMatrix[K][v] = mesh.connectivityMatrix[K][v];

    // Then handle nodes along edges
    for (int e = 0; e < mesh.numEdges; ++e)
    {
      const int& I = mesh.edgeArray[e][0];
      const int& J = mesh.edgeArray[e][1];

      // Grab nodes that form edge
      const MeshNode2D& A1 = mesh.meshNodes[I];
      const MeshNode2D& A2 = mesh.meshNodes[J];

      // Check if nodes form boundary edge
      bool boundary = false;
      if (mesh.edgeTypeMatrix[I][J] == 1)
        boundary = true;

      // Determine correct boundary condition
      BC_Type BC = BC_Type::Interior;
      if (A1.isCorner)
        BC = A2.BC;
      else
        BC = A1.BC;

      // Create p-1 equally spaced nodes along the line passing through A1,A2
      for (int i = 0; i < p - 1; ++i)
      {
        const real t = (i + 1.0) / p;
        const int n = mesh.numNodes + e * (p - 1) + i;
#pragma warning(suppress: 6386)
        FENodes[n].x = (1.0 - t) * A1.x + t * A2.x;
        FENodes[n].y = (1.0 - t) * A1.y + t * A2.y;
        if (boundary)
          FENodes[n].BC = BC;
#pragma warning(suppress: 6385)
        if ((int)FENodes[n].BC >= 0)
          boundaryIndices.emplace_back(n);
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
          const real tx = (real)(i + 1.0) / p;
          const real ty = (real)(j + 1.0) / p;

          // Transform points from reference domain to element K
          FENodes[nodeIndex].x = B[0][0] * tx + B[0][1] * ty + A1.x;
          FENodes[nodeIndex].y = B[1][0] * tx + B[1][1] * ty + A1.y;

          connectivityMatrix[K][localNodeIndex] = nodeIndex;
          ++nodeIndex;
          ++localNodeIndex;
        }
    }

    // Order boundary indices
    bool sorted = false;
    while (!sorted)
    {
      sorted = true;
      for (int i = 0; i < (int)boundaryIndices.size() - 1; ++i)
        if (boundaryIndices[i + 1] < boundaryIndices[i])
        {
          sorted = false;
          const int tmp = boundaryIndices[i];
          boundaryIndices[i] = boundaryIndices[i + 1];
          boundaryIndices[i + 1] = tmp;
        }
    }
  }

  /*
    Generates a finite element method using a given mesh, a
    polynomial order, and a (continuous) initial condition function.

    Will use Lagrange interpolation to approximate the initial condition function.
  */
  FEM2D(const Mesh2D& FEmesh, const int order, real2DFunction initialCondition)
    : FEM2D(FEmesh, order)
  {
    const auto& f = initialCondition;

    // Interpolate initialCondition
    for (int i = 0; i < Ng; ++i)
      FENodes[i][0] = f(FENodes[i].x, FENodes[i].y);
  }

  /*
    Move constructor.
  */
  FEM2D(FEM2D&& other) noexcept
    : mesh(other.mesh),
    polynomialOrder(other.polynomialOrder),
    Ng(other.Ng),
    boundaryIndices(std::move(other.boundaryIndices)),
    connectivityMatrix(std::move(connectivityMatrix))
  {
    FENodes = other.FENodes;
    other.FENodes = nullptr;
  }

  FEM2D(const FEM2D& other) = delete;

  FEM2D& operator=(const FEM2D& other) = delete;

  const Container<int>& operator[](const int elementIndex) const
  {
    // Debug
    ASSERT(elementIndex >= 0, "Element index must be non-negative");
    ASSERT(elementIndex < mesh.size, "Element index must be less than the number of elements");

    return connectivityMatrix[elementIndex];
  }

  /*
    \returns the FE node at the specified indices.

    \param nodeIndex: Index of the nodes within a particular element.
  */
  FENode2D<N>& operator()(const int elementIndex, const int nodeIndex) const
  {
    // Debug
    ASSERT(elementIndex >= 0, "Element index must be non-negative");
    ASSERT(elementIndex < mesh.size, "Element index must be less than the number of elements");
    ASSERT(nodeIndex >= 0, "Node index must be non-negative");
    ASSERT(nodeIndex < (polynomialOrder + 1)* (polynomialOrder + 2) / 2, "Node index must be less than the number of nodes per element");

    return FENodes[connectivityMatrix[elementIndex][nodeIndex]];
  }

  /*
    \returns the numerical approximation u_h(x, y).

    Note: This function is ineffiecient, as it searches for the element
    that (x, y) belongs to.  If the element is known, use other evaluate function.
  */
  real evaluate(const int varIndex, const real x, const real y) const { return evaluate(varIndex, x, y, 0, 0); }
  real evaluate(const int varIndex, const real x, const real y, const int xDerivativeOrder, const int yDerivativeOrder) const
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
    return evaluate(varIndex, x, y, K, xDerivativeOrder, yDerivativeOrder);
  }

  /*
    \returns the numerical approximation u_h(x, y) for a given
    derivative order.

    \param elementIndex: Index of element that contains (x, y).
  */
  real evaluate(const int varIndex, const real x, const real y, const int elementIndex, const int xDerivativeOrder, const int yDerivativeOrder) const
  {
    const int& K = elementIndex;
    const int& p = polynomialOrder;

    // Debug
    ASSERT(varIndex >= 0, "Variable index must be non-negative");
    ASSERT(varIndex < N, "Variables index must be less than the number of variables");
    ASSERT(K >= 0, "Element index must be non-negative");
    ASSERT(K < mesh.size, "Element index must be less than the number of elements");
    ASSERT(isInTriangle(x, y, K), "x must be in the element at the specified elementIndex");

    real sum = 0.0;
    for (int j = 0; j < (p + 1) * (p + 2) / 2; ++j)
    {
      const int& i = connectivityMatrix[K][j];
      sum += FENodes[i][varIndex] * lagrangeShapeFunction2D(x, y, *this, K, j, xDerivativeOrder, yDerivativeOrder);
    }
    return sum;
  }

  /*
    Outputs a file called "Plot.txt" that when graphed gives the
    d-th derivative of the FE interpolation.

    \param n: Number of points to plot per direction per element.
  */
  void plot(const int varIndex, const int n = 1) const { plot(varIndex, n, 0, 0); }
  void plot(const int varIndex, const int n, const int xDerivativeOrder, const int yDerivativeOrder) const
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

          file << x << ", " << y << ", " << evaluate(varIndex, x, y, K, xDerivativeOrder, yDerivativeOrder) << std::endl;
        }
    }
    file.close();
    std::cout << "Data written to \"Plot.txt\".  Press any key to continue" << std::endl;
    std::cin.get();
    remove("Plot.txt");
  }

  void plotField(const int n = 1)
  {
    ASSERT(N == 2, "Field plotting only supporting for two variable FEM structres");
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

          file << x << ", " << y << ", ";
          file << evaluate(0, x, y, K, 0, 0) << ", ";
          file << evaluate(1, x, y, K, 0, 0) << std::endl;
        }
    }
    file.close();
    std::cout << "Data written to \"Plot.txt\".  Press any key to continue" << std::endl;
    std::cin.get();
    remove("Plot.txt");
  }

  /*
    Determines if the given point (x, y) is inside the specified element.
  */
  bool isInTriangle(const real x, const real y, const int elementIndex) const
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

private:
  /*
    A 2D array that stores which FE nodes belong to each mesh element.
    connectivityMatrix[i] gives an array containing
    the indices of the nodes that belong to the i-th element.
  */
  Array2D<int> connectivityMatrix;
};