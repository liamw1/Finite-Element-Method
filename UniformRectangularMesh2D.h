#pragma once
#include "Precompilied.h"
#include "Mesh2D.h"
#include "Array2D.h"
#include "Print.h"

/*
  A rectangular 2D mesh consisting of uniform triangular elements.
*/
class UniformRectangularMesh2D : public Mesh2D
{
public:
  UniformRectangularMesh2D() = delete;

  /*
    Creates a 2(nx*ny) element mesh over the rectangular region
    xMin < x < xMax, yMin < y < yMax.

    \param nx: Number of elements along the x direction.
    \param ny: Number of elements along the y direction.
  */
  UniformRectangularMesh2D(const real xMin, const real xMax, const real yMin, const real yMax, const int nx, const int ny);

  UniformRectangularMesh2D(const UniformRectangularMesh2D& other) = delete;

  UniformRectangularMesh2D(UniformRectangularMesh2D&& other) noexcept;

  UniformRectangularMesh2D& operator=(UniformRectangularMesh2D& other) = delete;

  MeshNode2D operator()(const int elementIndex, const int nodeIndex) const;

private:
  const int numNodes;
  const int numEdges;

  /*
    A 2D array that stores which nodes belong to which mesh element.
    connectivityMatrix[i] gives an array containing
    the indices of the three nodes that belong to the i-th element.
  */
  Array2D<int> connectivityMatrix;

  /*
    A 2D array that stores which nodes belong to which mesh edge.
    edgeArray[i] gives an array containing
    the indices of the two nodes that form the i-th edge.
  */
  Array2D<int> edgeArray;

  /*
    A 2D array that stores information about the connections between nodes.

    edgeTypeMatrix[i][i] = p > 0 means the i-th node is a vertex of p mesh elements.
    edgeTypeMatrix[i][j] = 0 for i != j means nodes i and j do not form an edge.
    edgeTypeMatrix[i][j] = 1 for i != j means nodes i and j form a boundary edge.
    edgeTypeMatrix[i][j] = 2 for i != j means nodes i and j form an interior edge.

    Note: Full matrix representation is inefficient here.
    Better to use sparse matrix representation.
  */
  Array2D<int> edgeTypeMatrix;

  /*
    A 2D array that stores which nodes form which edges.

    edgeMatrix[i][j] = 0 means nodes i and j do not form an edge.
    edgeMatrix[i][j] = p > 0 means nodes i and j form the (p-1)-th edge.

    Note: Full matrix representation is inefficient here.
    Better to use sparse matrix representation.
  */
  Array2D<int> edgeMatrix;

  // Stores the x,y-values and boundary conditions of all the nodes in the mesh
  MeshNode2D* meshNodes;
};