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

template<int N>
struct FENode2D
{
  real x = 0.0;
  real y = 0.0;
  real var[N];
  BC_Type BC = BC_Type::Interior;
  bool isCorner = false;

  FENode2D()
  {
    ASSERT(N > 0, "Number of variables must be positive");
    for (int i = 0; i < N; ++i)
      var[i] = 0.0;
  }

  real& operator[](const int index)
  {
    ASSERT(index >= 0, "Variable index must be positive");
    ASSERT(index < N, "Variable index must be less than the number of variables");
    return var[index];
  }

  const real& operator[](const int index) const
  {
    ASSERT(index >= 0, "Variable index must be positive");
    ASSERT(index < N, "Variable index must be less than the number of variables");
    return var[index];
  }
};