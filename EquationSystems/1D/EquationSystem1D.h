#pragma once
#include "Precompilied.h"
#include "LinearAlgebra/Vector.h"
#include "LinearAlgebra/Matrix.h"
#include "Meshing/1D/FEM1D.h"

/*
  Interface for a general 1D boundary value problem equation system.
*/
class EquationSystem1D
{
public:
  virtual Vector solveSystem(const int n_gq) const = 0;

  virtual void update(const int n_gq) = 0;
};

Vector constructNaturalBoundaryVector1D(const FEM1D& fem, real1DFunction naturalBC);

Vector constructEssentialBoundaryVector1D(const FEM1D& fem, const Matrix& massMatrix);