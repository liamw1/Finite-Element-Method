#include "FEM2D.h"

FEM2D::FEM2D(const Mesh2D& FEmesh, const int order)
  : mesh(FEmesh), polynomialOrder(order)
{
}

FEM2D::FEM2D(FEM2D&& other) noexcept
  : mesh(other.mesh), polynomialOrder(other.polynomialOrder)
{
}
