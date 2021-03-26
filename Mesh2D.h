#pragma once
#include "Precompilied.h"
#include "BoundayEnums.h"
#include "Nodes.h"

/*
  Interface for a general 2D mesh.
*/
class Mesh2D
{
public:
  const int size;

  Mesh2D() = delete;

  Mesh2D(int n);

  Mesh2D(const Mesh2D& other) = delete;

  Mesh2D(Mesh2D&& other) noexcept;

  Mesh2D& operator=(const Mesh2D& other) = delete;

  virtual MeshNode2D operator()(const int elementIndex, const int nodeIndex) const = 0;

  virtual ~Mesh2D();
};