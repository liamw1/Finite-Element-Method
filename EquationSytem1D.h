#pragma once
#include "Precompilied.h"
#include "Vector.h"
#include "Matrix.h"
#include "FEM1D.h"

/*
  Interface for a general 1D boundary value problem equation system.
*/
class EquationSystem1D
{
public:
  virtual Vector solveSystem(const FEM1D& fem, const int n_gq) const = 0;
};

Vector constructNaturalBoundaryVector(const FEM1D& fem, real1DFunction naturalBC);

Vector constructEssentialBoundaryVector(const FEM1D& fem, const Matrix& massMatrix);