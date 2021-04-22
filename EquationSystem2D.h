#pragma once
#include "Precompilied.h"
#include "Vector.h"
#include "Matrix.h"
#include "FEM2D.h"
#include "Gauss-LegendreNodes.h"
#include "LagrangeShapeFunctions2D.h"

/*
  Interface for a general 2D boundary value problem equation system.
*/
class EquationSystem2D
{
public:
  virtual Vector solveSystem(const FEM2D& fem, const int n_gq) const = 0;

protected:
  void removeBoundaryIndices(Vector& v, const std::vector<int> bI) const;
  void removeBoundaryIndices(Matrix& A, const std::vector<int> bI) const;
};

Vector constructEssentialBoundaryVector2D(const FEM2D& fem, const Matrix& massMatrix);

Vector constructNaturalBoundaryVector2D(const FEM2D& fem, real2DFunction naturalBC, const int n_gq);