#include "Precompilied.h"
#include "EquationSystem2D.h"

Vector constructEssentialBoundaryVector2D(const FEM2D& fem, const Matrix& massMatrix)
{
  const int& p = fem.polynomialOrder;

  Vector bc_e = Vector(fem.Ng);
  for (int i = 0; i < fem.Ng; ++i)
    for (int n = 0; n < fem.boundaryIndices.size(); ++n)
    {
      const int& j = fem.boundaryIndices[n];
      bc_e[i] -= massMatrix[i][j] * fem.FENodes[j].u;
    }
  return bc_e;
}

Vector constructNaturalBoundaryVector2D(const FEM2D& fem, real2DFunction naturalBC, const int n_gq)
{
  const int& p = fem.polynomialOrder;

  Vector bc_n = Vector(fem.Ng);
  for (int K = 0; K < fem.mesh.size; ++K)
  {
    // Find boundary edge if element has one
    int e = -1;
    for (int i = 0; i < 3; ++i)
      for (int j = i + 1; j < 3; ++j)
      {
        const int I = fem.mesh.connectivityMatrix[K][i];
        const int J = fem.mesh.connectivityMatrix[K][j];
        if (fem.mesh.edgeTypeMatrix[I][J] == 1)
        {
          ASSERT(fem.mesh.edgeMatrix[I][J] > 0, "Nodes do not form an edge");
          e = fem.mesh.edgeMatrix[I][J] - 1;
        }
      }

    if (e >= 0)
    {
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
          bc_n[fem[K][j]] += innerProduct; // Accumulate to bc_n
        }
      }
    }
  }
  return bc_n;
}

void EquationSystem2D::removeBoundaryIndices(Vector& v, const std::vector<int> bI) const
{
  if (bI.size() > 0)
  {
    Vector u = Vector(v.size() - bI.size());
    int index = 0;
    for (int i = 0; i < v.size(); ++i)
    {
      if (i != bI[index])
        u[i - index] = v[i];
      else
        ++index;
    }
    v = std::move(u);
  }
}

void EquationSystem2D::removeBoundaryIndices(Matrix& A, const std::vector<int> bI) const
{
  if (bI.size() > 0)
  {
    Matrix M = Matrix(A.size() - bI.size());
    int rowIndex = 0;
    for (int i = 0; i < A.size(); ++i)
    {
      if (i != bI[rowIndex])
      {
        int columnIndex = 0;
        for (int j = 0; j < A.size(); ++j)
        {
          if (j != bI[columnIndex])
            M[i - rowIndex][j - columnIndex] = A[i][j];
          else
            ++columnIndex;
        }
      }
      else
        ++rowIndex;
    }
    A = std::move(M);
  }
}
