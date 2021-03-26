#include "Mesh2D.h"

Mesh2D::Mesh2D(int n)
  : size(n)
{
}

Mesh2D::Mesh2D(Mesh2D&& other) noexcept
  : size(other.size)
{
}

Mesh2D::~Mesh2D()
{
}
