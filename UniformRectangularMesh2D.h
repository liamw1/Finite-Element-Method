#pragma once
#include "Precompilied.h"
#include "Mesh2D.h"
#include "Array2D.h"
#include "Print.h"

/*
  A rectangular 2D mesh consisting of uniform triangular elements.
*/
class UniformRectangularMesh2D : public Mesh2D
{
public:
  UniformRectangularMesh2D() = delete;

  /*
    Creates a 2(nx*ny) element mesh over the rectangular region
    xMin < x < xMax, yMin < y < yMax.

    \param nx: Number of elements along the x direction.
    \param ny: Number of elements along the y direction.
  */
  UniformRectangularMesh2D(const real xMin, const real xMax, const real yMin, const real yMax, const int nx, const int ny);

  UniformRectangularMesh2D(const UniformRectangularMesh2D& other) = delete;

  UniformRectangularMesh2D(UniformRectangularMesh2D&& other) noexcept;

  UniformRectangularMesh2D& operator=(UniformRectangularMesh2D& other) = delete;

  MeshNode2D operator()(const int elementIndex, const int nodeIndex) const;

private:
};