#pragma once
#include "Precompilied.h"
#include "Mesh2D.h"
#include "Array2D.h"

/*
  Generates a general 2D mesh from a file.

  File should have two sections: 

    1) A list of node coordinates beginning with the 
    line ">StartNodes" and ending with ">EndNodes".

    2) A list of which nodes belong to which elements
    (to form connectivity matrix), beginning with the
    line ">StartElements" and ending with ">EndElements".

  Note: Boundary conditions have not been implemented for this mesh!
*/
class UnstructuredMesh2D : public Mesh2D
{
public:
  UnstructuredMesh2D() = delete;

  /*
    Constructs mesh from file with specified name.
  */
  UnstructuredMesh2D(const std::string meshFile);

  UnstructuredMesh2D(const UnstructuredMesh2D& other) = delete;

  UnstructuredMesh2D(UnstructuredMesh2D&& other) noexcept;

  UnstructuredMesh2D& operator=(const UnstructuredMesh2D& other) = delete;

  MeshNode2D operator()(const int elementIndex, const int nodeIndex) const override;

private:
};