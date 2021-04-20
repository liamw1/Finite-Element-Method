#pragma once
#include "Precompilied.h"
#include "EquationSystem1D.h"
#include "EquationSystem2D.h"

/*
  Updates solution based on equation system provided.

  Note: Essential boundary conditions should be enforeced
  before each call of this function.
*/
void update1D(FEM1D& fem, const EquationSystem1D& equationSystem, const int n_gq);

/*
  Updates solution based on equation system provided.

  Note: Essential boundary conditions should be enforeced
  before each call of this function.
*/
void update2D(FEM2D& fem, const EquationSystem2D& equationSystem, const int n_gq);