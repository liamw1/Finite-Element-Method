#pragma once
#include "Precompilied.h"
#include "BoundayEnums.h"

struct MeshNode
{
  real x = 0.0;
  BC_Type BC = BC_Type::Interior;
};

struct FENode
{
  real x = 0.0;
  real u = 0.0;
  BC_Type BC = BC_Type::Interior;
};