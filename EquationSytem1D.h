#pragma once
#include "Precompilied.h"
#include "Vector.h"
#include "Matrix.h"
#include "FEM.h"

/*
  Interface for a general 1D boundary value problem equation system.
*/
class EquationSystem1D
{
public:
  virtual Vector solveSystem(const FEM& fem, const int n_gq) const = 0;
};

Vector constructNaturalBoundaryVector(const FEM& fem, realFunction naturalBC);

Vector constructEssentialBoundaryVector(const FEM& fem, const Matrix& massMatrix);