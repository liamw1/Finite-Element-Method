#pragma once
#include "Precompilied.h"
#include "BoundayEnums.h"
#include "Nodes.h"

/*
  Interface for a general 1D mesh.
*/
class Mesh1D
{
public:
  const int size;  // Number of mesh elements
  const real xL, xR;
  int numBoundaryNodes = -1;

  Mesh1D() = delete;

  Mesh1D(const real leftEndPoint, const real rightEndPoint, const int n);

  Mesh1D(const Mesh1D& other) = delete;

  Mesh1D(Mesh1D&& other) noexcept;

  Mesh1D& operator=(const Mesh1D& other) = delete;

  virtual MeshNode1D operator()(const int elementIndex, const EdgeType edgeIndex) const = 0;

  virtual void setBoundaryConditions(const BC_Type leftBoundaryCondition, const BC_Type rightBoundaryCondition) = 0;

  virtual ~Mesh1D();
};