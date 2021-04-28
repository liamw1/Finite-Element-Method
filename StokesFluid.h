#pragma once
#include "Precompilied.h"
#include "EquationSystem2D.h"
#include "L2Projection.h"

/*
  Equation class for a steady state Stokes fluid BVP.
*/
class StokesFluid : public EquationSystem2D
{
public:
  StokesFluid() = delete;

  StokesFluid(FEM2D<2>& uFEM, FEM2D<1>& pFEM, real2DFunction fFunc, real2DFunction nuFunc, real2DFunction rhoFunc, real2DFunction naturalBoundaryCondition);

  const int neq() const override;

  /*
    \returns coefficients of FE approximation u_h.

    \param n_gq: Number of Gaussian quadrature nodes.

    Essential boundary conditions should be enforeced
    before each call of this function.
  */
  Vector solveSystem(const int n_gq) const override;

private:
  real2DFunction f;          // Body force
  real2DFunction nu;         // Kinematic viscosity
  real2DFunction rho;        // Density
  real2DFunction naturalBC;

  FEM2D<2>& uFem;
  FEM2D<1>& pFem;
};