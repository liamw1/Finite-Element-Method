#include "Precompilied.h"
#include "Hwk4_C1.h"
#include "Hwk4_C2.h"
#include "Hwk4_C3.h"
#include "Hwk6_B2.h"
#include "Hwk8_C1.h"
#include "L2Projection.h"
#include "UniformRectangularMesh2D.h"
#include "UnstructuredMesh2D.h"

std::ostream& operator<<(std::ostream& os, const BC_Type BC)
{
  os << (int)BC;
  return os;
}

int main()
{
  UniformRectangularMesh2D mesh = UniformRectangularMesh2D(0, 1, 0, 1, 5, 5);
  mesh.setBoundaryConditions(BC_Type::Natural, BC_Type::Natural, BC_Type::Natural, BC_Type::Natural);
  for (int i = 0; i < mesh.numNodes; ++i)
  {
    print("(", mesh.meshNodes[i].x, ", ", mesh.meshNodes[i].y, ")");
    print(mesh.meshNodes[i].BC);
    print();
  }
}