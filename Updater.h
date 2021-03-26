#pragma once
#include "Precompilied.h"
#include "EquationSytem1D.h"

/*
  Updates solution based on equation system provided.

  Essential boundary conditions should be enforeced
  before each call of this function.
*/
void update(FEM1D& fem, const EquationSystem1D& equationSystem, const int n_gq);