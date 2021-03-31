#pragma once
#include "Precompilied.h"
#include "Mesh2D.h"
#include "Matrix.h"

/*
  2D finite element structure.

  Assumes triangular mesh.
*/
class FEM2D
{
public:
  const Mesh2D& mesh;
  const int polynomialOrder;
  const int Ng;
  FENode2D* FENodes;

  FEM2D() = delete;

  /*
    Generates a finite element method using a given mesh and a
    polynomial order.
  */
  FEM2D(const Mesh2D& FEmesh, const int order);

  /*
    Move constructor.
  */
  FEM2D(FEM2D&& other) noexcept;

  FEM2D(const FEM2D& other) = delete;

  FEM2D& operator=(const FEM2D& other) = delete;

private:
  /*
    A 2D array that stores which FE nodes belong to each mesh element.
    connectivityMatrix[i] gives an array containing
    the indices of the nodes that belong to the i-th element.
  */
  Array2D<int> connectivityMatrix;
};