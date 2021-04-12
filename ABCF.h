#pragma once
#include "Precompilied.h"
#include "EquationSytem1D.h"
#include "L2Projection.h"

/*
  Equation class for BVP of the form -(au')' + bu' + cu = f.  With boundary conditions
  a(x1)u'(x1) = g0 and u(x2) = g1, where x1, x2 are endpoints of the domain.
*/
class ABCF : public EquationSystem1D
{
public:
  ABCF() = delete;

  ABCF(real1DFunction aFunc, real1DFunction bFunc, real1DFunction cFunc, real1DFunction fFunc, real1DFunction naturalBoundaryCondition);

  /*
    \returns coefficients of FE approximation u_h.

    \param n_gq: Number of Gaussian quadrature nodes.

    Essential boundary conditions should be enforeced
    before each call of this function.
  */
  Vector solveSystem(const FEM1D& fem, const int n_gq) const override;

private:
  real1DFunction a, b, c, f;
  real1DFunction naturalBC;
};