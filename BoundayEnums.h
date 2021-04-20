#pragma once
#include "Precompilied.h"

enum class EdgeType
{
  Left,
  Right
};

/*
  Flag for FE nodes based on the boundary condition that will be imposed on them.

  Nodes marked with "Interior" should not be on the boundary 
  and should not have any boundary conditions.
*/
enum class BC_Type
{
  Interior = -200,  // Default value
  Natural = -1,
  Essential = 0,
  Dirichlet = 1,
  Neumann = 2,
  Corner = 200
};

enum class MeshType
{
  Triangular = 3
};