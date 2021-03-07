#include "Precompilied.h"
#include "UniformMesh1D.h"

UniformMesh1D::UniformMesh1D(const real leftEndPoint, const real rightEndPoint, int n)
  : Mesh1D(leftEndPoint, rightEndPoint, n), connectivityMatrix(Array2D<int>(n, 2))
{
  // Debug
  ASSERT(n > 0, "Invalid mesh size: A mesh must have a least one element!");
  ASSERT(xL < xR, "Invalid domain: xR must be greater than xL!");

  meshNodes = new MeshNode[n + 1];
  meshNodes[0].x = xL;
  for (int i = 0; i < n; ++i)
  {
#pragma warning(suppress: 6386)
    meshNodes[i + 1].x = xL + (i + 1.0) * (xR - xL) / n;
    connectivityMatrix[i][(int)EdgeType::Left] = i;
    connectivityMatrix[i][(int)EdgeType::Right] = i + 1;
  }
}

UniformMesh1D::UniformMesh1D(UniformMesh1D&& other) noexcept
  : Mesh1D(other.xL, other.xR, other.size), connectivityMatrix(std::move(other.connectivityMatrix))
{
  meshNodes = other.meshNodes;
  other.meshNodes = nullptr;
}

MeshNode UniformMesh1D::operator()(const int elementIndex, const EdgeType edgeIndex) const
{
  // Debug
  ASSERT(elementIndex >= 0, "Element index must be non-negative");
  ASSERT(elementIndex < size, "Element index must be less than the number of elements");

  return meshNodes[connectivityMatrix[elementIndex][(int)edgeIndex]];
}

UniformMesh1D::~UniformMesh1D()
{
  delete[] meshNodes;
}

void UniformMesh1D::setBoundaryConditions(const BC_Type leftBoundaryCondition, const BC_Type rightBoundaryCondition)
{
  meshNodes[0].BC = leftBoundaryCondition;
  meshNodes[size].BC = rightBoundaryCondition;

  if ((int)leftBoundaryCondition >= 0 && (int)rightBoundaryCondition >= 0)
    numBoundaryNodes = 2;
  else if ((int)leftBoundaryCondition >= 0 || (int)rightBoundaryCondition >= 0)
    numBoundaryNodes = 1;
}
