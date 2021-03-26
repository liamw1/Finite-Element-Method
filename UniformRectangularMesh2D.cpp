#include "Precompilied.h"
#include "UniformRectangularMesh2D.h"

UniformRectangularMesh2D::UniformRectangularMesh2D(const real xMin, const real xMax, const real yMin, const real yMax, const int nx, const int ny)
  : Mesh2D(2 * nx * ny), connectivityMatrix(Array2D<int>(size, 3))
{
  // Debug
  ASSERT(size > 0, "Invalid mesh size: A mesh must have a least one element!");

  const real dx = (xMax - xMin) / nx;
  const real dy = (yMax - yMin) / ny;

  // Set coordinates of nodes
  meshNodes = new MeshNode2D[(nx + 1) * (ny + 1)];
  for (int i = 0; i < ny + 1; ++i)
    for (int j = 0; j < nx + 1; ++j)
    {
#pragma warning(suppress: 6386)
      meshNodes[i * (nx+1) + j].x = xMin + j * dx;
      meshNodes[i * (nx+1) + j].y = yMin + i * dy;
    }

  int K = 0;
  // Loop over rectangular regions
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
      ++K;

      connectivityMatrix[K][0] = botLeft;
      connectivityMatrix[K][1] = topLeft;
      connectivityMatrix[K][2] = topRight;
      ++K;
    }
}

UniformRectangularMesh2D::UniformRectangularMesh2D(UniformRectangularMesh2D&& other) noexcept
  : Mesh2D(size), connectivityMatrix(std::move(other.connectivityMatrix))
{
  meshNodes = other.meshNodes;
  other.meshNodes = nullptr;
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
