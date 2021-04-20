#include "Precompilied.h"
#include "UniformRectangularMesh2D.h"

UniformRectangularMesh2D::UniformRectangularMesh2D(const real xMin, const real xMax, const real yMin, const real yMax, const int nx, const int ny)
  : Mesh2D(2 * nx * ny, (nx + 1)* (ny + 1), 5 + 4 * (nx + ny - 2) + 3 * (nx - 1) * (ny - 1)),
    xL(xMin), xR(xMax), yL(yMin), yR(yMax)
{
  // Debug
  ASSERT(size > 0, "Invalid mesh size: A mesh must have a least one element!");

  // Initialize data structures
  meshNodes = new MeshNode2D[numNodes];
  connectivityMatrix = Array2D<int>(size, 3);
  edgeArray = Array2D<int>(numEdges, 2);
  edgeTypeMatrix = Array2D<int>(numNodes, numNodes);
  edgeMatrix = Array2D<int>(numNodes, numNodes);

  const real dx = (xMax - xMin) / nx;
  const real dy = (yMax - yMin) / ny;

  // Set node coordinates
  for (int i = 0; i < ny + 1; ++i)
    for (int j = 0; j < nx + 1; ++j)
    {
#pragma warning(suppress: 6386)
      meshNodes[i * (nx+1) + j].x = xMin + j * dx;
      meshNodes[i * (nx+1) + j].y = yMin + i * dy;
    }
  
  // Form connectivity matrix
  int K = 0;
  // Slice region into dx*dy rectangles and loop over them
  for (int i = 0; i < ny; ++i)
    for (int j = 0; j < nx; ++j)
    {
      // Indices corresponding to each corner
      const int botLeft = i * (nx + 1) + j;
      const int botRight = i * (nx + 1) + j + 1;
      const int topLeft = (i + 1) * (nx + 1) + j;
      const int topRight = (i + 1) * (nx + 1) + j + 1;

      // Bisect rectangle to form two triangular mesh elements
      connectivityMatrix[K][0] = botLeft;
      connectivityMatrix[K][1] = botRight;
      connectivityMatrix[K][2] = topRight;
      connectivityMatrix[K + 1][0] = botLeft;
      connectivityMatrix[K + 1][1] = topLeft;
      connectivityMatrix[K + 1][2] = topRight;
      K += 2;
    }

  // Form edge type matrix
  for (int i = 0; i < numNodes; ++i)
    for (int j = 0; j < numNodes; ++j)
      edgeTypeMatrix[i][j] = 0;
  for (int i = 0; i < size; ++i)
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 3; ++k)
      {
        // Grabbing indices of j-th and k-th vertex in element
        const int& vj = connectivityMatrix[i][j];
        const int& vk = connectivityMatrix[i][k];

        ++edgeTypeMatrix[vj][vk];
      }

  // Form edge array
  int E = 0;
  for (int i = 0; i < numNodes; ++i)
    for (int j = i + 1; j < numNodes; ++j)
    {
      if (edgeTypeMatrix[i][j] > 0)
      {
        edgeArray[E][0] = i;
        edgeArray[E][1] = j;
        ++E;
      }
    }

  // Form edge matrix
  for (int i = 0; i < numNodes; ++i)
    for (int j = 0; j < numNodes; ++j)
      edgeMatrix[i][j] = 0;
  for (int i = 0; i < numEdges; ++i)
  {
    edgeMatrix[edgeArray[i][0]][edgeArray[i][1]] = i + 1;
    edgeMatrix[edgeArray[i][1]][edgeArray[i][0]] = i + 1;
  }
}

UniformRectangularMesh2D::UniformRectangularMesh2D(UniformRectangularMesh2D&& other) noexcept
  : Mesh2D(other.size, other.numNodes, other.numEdges),
    xL(other.xL), xR(other.xR), yL(other.yL), yR(other.yR)
{
}

MeshNode2D UniformRectangularMesh2D::operator()(const int elementIndex, const int nodeIndex) const
{
  // Debug
  ASSERT(elementIndex >= 0, "Element index must be non-negative");
  ASSERT(elementIndex < size, "Element index must be less than the number of elements");
  ASSERT(nodeIndex >= 0, "Node index must be non-negative");
  ASSERT(nodeIndex < 3, "Node index must be less than 3 on a triangular mesh");

  return meshNodes[connectivityMatrix[elementIndex][nodeIndex]];
}

void UniformRectangularMesh2D::setBoundaryConditions(const BC_Type boundaryCondition)
{
  setBoundaryConditions(boundaryCondition, boundaryCondition, boundaryCondition, boundaryCondition);
}

void UniformRectangularMesh2D::setBoundaryConditions(const BC_Type leftBoundaryCondition,
                                                     const BC_Type rightBoundaryCondition,
                                                     const BC_Type bottomBoundaryCondition,
                                                     const BC_Type topBoundaryCondition)
{
  // Find edge nodes
  for (int i = 0; i < numNodes; ++i)
  {
    if (meshNodes[i].x == xL)
      meshNodes[i].BC = leftBoundaryCondition;
    if (meshNodes[i].x == xR)
      meshNodes[i].BC = rightBoundaryCondition;
    if (meshNodes[i].y == yL)
      meshNodes[i].BC = bottomBoundaryCondition;
    if (meshNodes[i].y == yR)
      meshNodes[i].BC = topBoundaryCondition;
  }

  // Find corner nodes
  for (int i = 0; i < numNodes; ++i)
    if ((meshNodes[i].x == xL || meshNodes[i].x == xR) && (meshNodes[i].y == yL || meshNodes[i].y == yR))
      meshNodes[i].BC = BC_Type::Corner;

  // Count number of boundary nodes
  numBoundaryNodes = 0;
  for (int i = 0; i < numNodes; ++i)
    if ((int)meshNodes[i].BC >= 0)
      ++numBoundaryNodes;
}
