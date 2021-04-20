#pragma once
#include "Precompilied.h"
#include "BoundayEnums.h"
#include "Nodes.h"
#include "Array2D.h"

/*
  Interface for a general 2D mesh.
*/
class Mesh2D
{
public:
  const int size;           // Number of mesh elements
  const int numNodes;       // Number of mesh nodes
  const int numEdges;       // Number of mesh edges
  int numBoundaryNodes = -1;     // Number of boundary nodes

  /*
    A 2D array that stores which nodes belong to which mesh element.
    connectivityMatrix[i] gives an array containing
    the indices of the three nodes that belong to the i-th element.
  */
  Array2D<int> connectivityMatrix{};

  /*
    A 2D array that stores which nodes belong to which mesh edge.
    edgeArray[i] gives an array containing
    the indices of the two nodes that form the i-th edge.
  */
  Array2D<int> edgeArray{};

  /*
    A 2D array that stores information about the connections between nodes.

    edgeTypeMatrix[i][i] = p > 0 means the i-th node is a vertex of p mesh elements.
    edgeTypeMatrix[i][j] = 0 for i != j means nodes i and j do not form an edge.
    edgeTypeMatrix[i][j] = 1 for i != j means nodes i and j form a boundary edge.
    edgeTypeMatrix[i][j] = 2 for i != j means nodes i and j form an interior edge.

    Note: Full matrix representation is inefficient here.
    Better to use sparse matrix representation.
  */
  Array2D<int> edgeTypeMatrix{};

  /*
    A 2D array that stores which nodes form which edges.

    edgeMatrix[i][j] = 0 means nodes i and j do not form an edge.
    edgeMatrix[i][j] = p > 0 means nodes i and j form the (p-1)-th edge.

    Note: Full matrix representation is inefficient here.
    Better to use sparse matrix representation.
  */
  Array2D<int> edgeMatrix{};

  // Stores the x,y-values and boundary conditions of all the nodes in the mesh
  MeshNode2D* meshNodes = nullptr;

  Mesh2D() = delete;

  Mesh2D(const int n, const int nodes, const int edges);

  Mesh2D(const Mesh2D& other) = delete;

  Mesh2D(Mesh2D&& other) noexcept;

  Mesh2D& operator=(const Mesh2D& other) = delete;

  virtual MeshNode2D operator()(const int elementIndex, const int nodeIndex) const = 0;

  virtual ~Mesh2D();
};