#pragma once
#include "Precompilied.h"
#include "Mesh1D.h"
#include "Array2D.h"

/*
  A 1D mesh consisting of n equally spaces nodes.
*/
class UniformMesh1D : public Mesh1D
{
public:

  UniformMesh1D() = delete;

  /*
    Generates a uniform 1D mesh with n elements on the interval [xL, xR].
  */
  UniformMesh1D(const real leftEndpoint, const real rightEndPoint, int n);

  UniformMesh1D(const UniformMesh1D& other) = delete;

  /*
    Move constructor.
  */
  UniformMesh1D(UniformMesh1D&& other) noexcept;

  UniformMesh1D& operator=(const UniformMesh1D& other) = delete;

  /*
    \returns The node at the specified element index
    and egdeIndex.

    \param edgeIndex: Set to 0 for left edge and 1 for right edge.
  */
  MeshNode operator()(const int elementIndex, const EdgeType edgeIndex) const override;

  ~UniformMesh1D() override;

  void setBoundaryConditions(const BC_Type leftBoundaryCondition, const BC_Type rightBoundaryCondition) override;

private:
  /*
    A 2D array that stores which edge nodes belong to each mesh element.
    connectivityMatrix[i] gives an array containing
    the indices of the two nodes that belong to the i-th element.
  */
  Array2D<int> connectivityMatrix;

  // Stores the x-values and boundary conditions of all the nodes in the mesh
  MeshNode* meshNodes;
};