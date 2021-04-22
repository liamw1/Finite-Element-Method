#pragma once
#include "Precompilied.h"
#include "BoundayEnums.h"

struct MeshNode1D
{
  real x = 0.0;
  BC_Type BC = BC_Type::Interior;
};

struct FENode1D
{
  real x = 0.0;
  real u = 0.0;
  BC_Type BC = BC_Type::Interior;
};

struct MeshNode2D
{
  real x = 0.0;
  real y = 0.0;
  BC_Type BC = BC_Type::Interior;
  bool isCorner = false;
};

struct FENode2D
{
  real x = 0.0;
  real y = 0.0;
  real u = 0.0;
  BC_Type BC = BC_Type::Interior;
  bool isCorner = false;
};