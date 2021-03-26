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

  UniformRectangularMesh2D(const real xMin, const real xMax, const real yMin, const real yMax, const int nx, const int ny);

  UniformRectangularMesh2D(const UniformRectangularMesh2D& other) = delete;

  UniformRectangularMesh2D(UniformRectangularMesh2D&& other) noexcept;

  UniformRectangularMesh2D& operator=(UniformRectangularMesh2D& other) = delete;

  MeshNode2D operator()(const int elementIndex, const int nodeIndex) const;

private:
  /*
    A 2D array that stores which nodes belong to which mesh element.
    connectivityMatrix[i] gives an array containing
    the indices of the three nodes that belong to the i-th element.
  */
  Array2D<int> connectivityMatrix;

  // Stores the x,y-values and boundary conditions of all the nodes in the mesh
  MeshNode2D* meshNodes;
};