#pragma once
#include "Precompilied.h"
#include "Mesh2D.h"

class FEM2D
{
public:
  const Mesh2D& mesh;
  const int polynomialOrder;

  FEM2D() = delete;

  /*
    Generates a finite element method using a given mesh and a
    polynomial order.
  */
  FEM2D(const Mesh2D& FEmesh, const int order);

  /*
    Move constructor.
  */
  FEM2D(FEM2D&& other) noexcept;

  FEM2D(const FEM2D& other) = delete;

  FEM2D& operator=(const FEM2D& other) = delete;

private:

};