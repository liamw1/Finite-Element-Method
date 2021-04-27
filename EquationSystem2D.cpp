#include "Precompilied.h"
#include "EquationSystem2D.h"

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
