#pragma once
#include "Precompilied.h"
#include "EquationSytem1D.h"
#include "L2Projection.h"

/*
  Equation class for BVP of the form -(au')' + cu = f.  With boundary conditions
  a(x1)u'(x1) = g0 and u(x2) = g1, where x1, x2 are endpoints of the domain.
*/
class ACF : public EquationSystem1D
{
public:
  ACF() = delete;

  ACF(realFunction aFunc, realFunction cFunc, realFunction fFunc, realFunction naturalBoundaryCondition);

  /*
    \returns coefficients of FE approximation u_h.

    \param n_gq: Number of Gaussian quadrature nodes.

    Essential boundary conditions should be enforeced
    before each call of this function.
  */
  Vector solveSystem(const FEM1D& fem, const int n_gq) const override;

private:
  realFunction a, c, f;
  realFunction naturalBC;
};