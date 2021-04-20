#pragma once
#include "Precompilied.h"
#include "EquationSystem2D.h"
#include "L2Projection.h"

/*
  Equation class for BVP of the form -div(a*grad(u)) + cu = f.  With the boundary
  condition u|dOmega_D = g_D.
*/
class Elliptic2DACF : public EquationSystem2D
{
public:
  Elliptic2DACF() = delete;

  Elliptic2DACF(real2DFunction aFunc, real2DFunction cFunc, real2DFunction fFunc);

  /*
    \returns coefficients of FE approximation u_h.

    \param n_gq: Number of Gaussian quadrature nodes.

    Essential boundary conditions should be enforeced
    before each call of this function.
  */
  Vector solveSystem(const FEM2D& fem, const int n_gq) const override;

private:
  real2DFunction a, c, f;
};