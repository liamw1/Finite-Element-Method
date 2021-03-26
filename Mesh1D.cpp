#include "Precompilied.h"
#include "Mesh1D.h"

Mesh1D::Mesh1D(const real leftEndPoint, const real rightEndPoint, const int n)
  : xL(leftEndPoint), xR(rightEndPoint), size(n)
{
}

Mesh1D::Mesh1D(Mesh1D&& other) noexcept
  : xL(other.xL), xR(other.xR), size(other.size)
{
}

Mesh1D::~Mesh1D()
{
}
