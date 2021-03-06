#pragma once
#include "Precompilied.h"
#include "LinearAlgebra/Vector.h"
#include "LinearAlgebra/Matrix.h"
#include "Meshing/2D/FEM2D.h"
#include "Functions/Gauss-LegendreNodes.h"
#include "Functions/LagrangeShapeFunctions2D.h"

/*
  Interface for a general 2D boundary value problem equation system.
*/
class EquationSystem2D
{
public:
  /*
    Number of equations/variables in the system
  */
  virtual int neq() const = 0;

  virtual Vector solveSystem(const int n_gq) const = 0;

  virtual void update(const int n_gq) = 0;

protected:
  void removeBoundaryIndices(Vector* v, const std::vector<int> bI) const;
  void removeBoundaryIndices(Matrix* A, const std::vector<int> bI) const;
  void removeBoundaryIndices(Matrix* A, const std::vector<int> rI, const std::vector<int> cI) const;
};

template<int N>
Vector* constructEssentialBoundaryVector2D(const FEM2D<N>& fem, const int varIndex, const Matrix* massMatrix)
{
  return constructEssentialBoundaryVector2D(fem, fem, varIndex, massMatrix);
}
template<int N, int M>
Vector* constructEssentialBoundaryVector2D(const FEM2D<N>& fem1, const FEM2D<M>& fem2, const int varIndex, const Matrix* massMatrix)
{
  Vector* bc_e = new Vector(fem1.Ng);
  for (int i = 0; i < fem1.Ng; ++i)
    for (int n = 0; n < fem2.boundaryIndices.size(); ++n)
    {
      const int& j = fem2.boundaryIndices[n];
      (*bc_e)[i] -= (*massMatrix)[i][j] * fem2.FENodes[j][varIndex];
    }
  return bc_e;
}

template<int N>
Vector* constructNaturalBoundaryVector2D(const FEM2D<N>& fem, real2DFunction naturalBC, const int n_gq)
{
  const int& p = fem.polynomialOrder;

  Vector* bc_n = new Vector(fem.Ng);
  for (int K = 0; K < fem.mesh.size; ++K)
  {
    // Find boundary edge if element has one
    for (int i = 0; i < 3; ++i)
      for (int j = i + 1; j < 3; ++j)
      {
        const int I = fem.mesh.connectivityMatrix[K][i];
        const int J = fem.mesh.connectivityMatrix[K][j];
        if (fem.mesh.edgeTypeMatrix[I][J] == 1)
        {
          ASSERT(fem.mesh.edgeMatrix[I][J] > 0, "Nodes do not form an edge");
          const int e = fem.mesh.edgeMatrix[I][J] - 1;

          // Check to see if boundary conditions should be applied
          const MeshNode2D& A1 = fem.mesh.meshNodes[fem.mesh.edgeArray[e][0]];
          const MeshNode2D& A2 = fem.mesh.meshNodes[fem.mesh.edgeArray[e][1]];
          ASSERT(A1.BC != BC_Type::Interior && A2.BC != BC_Type::Interior, "Nodes do not form boundary edge");
          if (A1.BC == BC_Type::Natural || A2.BC == BC_Type::Natural)
          {
            // Grab 1D quadrature nodes along edge
            std::vector<std::array<real, 2>> GLnodes = gaussEdgeNodesLocal(fem.mesh, e, n_gq);
            std::vector<real> GLweights = gaussEdgeWeightsLocal(fem.mesh, e, n_gq);

            for (int j = 0; j < (p + 1) * (p + 2) / 2; ++j)
            {
              // Calculate inner product between naturalBC and j-th shape function on K
              real innerProduct = 0.0;
              for (int i = 0; i < n_gq; ++i)
              {
                const real& x = GLnodes[i][0];
                const real& y = GLnodes[i][1];
                const real integrand = naturalBC(x, y) * lagrangeShapeFunction2D(x, y, fem, K, j, 0, 0);
                innerProduct += GLweights[i] * integrand;
              }
              (*bc_n)[fem[K][j]] += innerProduct; // Accumulate to bc_n
            }
          }
        }
      }
  }
  return bc_n;
}