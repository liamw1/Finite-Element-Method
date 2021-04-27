#pragma once
#include "Precompilied.h"
#include "EquationSystem2D.h"
#include "L2Projection.h"

/*
  Equation class for BVP of the form -div(a*grad(u)) + b*div(u) + cu = f.  With the boundary
  condition u|dOmega_D = g_D.
*/
class Elliptic2DABCF : public EquationSystem2D
{
public:
  Elliptic2DABCF() = delete;

  Elliptic2DABCF(FEM2D<1>& uFem, real2DFunction aFunc, real2DFunction bFunc, real2DFunction cFunc, real2DFunction fFunc, real2DFunction naturalBoundaryCondition);

  const int neq() const override;

  /*
    \returns coefficients of FE approximation u_h.

    \param n_gq: Number of Gaussian quadrature nodes.

    Essential boundary conditions should be enforeced
    before each call of this function.
  */
  Vector solveSystem(const int n_gq) const override;

private:
  FEM2D<1>& fem;
  real2DFunction a, b, c, f;
  real2DFunction naturalBC;
};