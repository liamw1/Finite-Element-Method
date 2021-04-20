#pragma once
#include "Precompilied.h"
#include "Vector.h"
#include "Matrix.h"
#include "FEM2D.h"

/*
  Interface for a general 2D boundary value problem equation system.
*/
class EquationSystem2D
{
public:
  virtual Vector solveSystem(const FEM2D& fem, const int n_gq) const = 0;
};

Vector constructEssentialBoundaryVector2D(const FEM2D& fem, const Matrix& massMatrix);