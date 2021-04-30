#include "Precompilied.h"
#include "EquationSystem2D.h"

void EquationSystem2D::removeBoundaryIndices(Vector& v, const std::vector<int> bI) const
{
  // Debug
  ASSERT(v.size() >= bI.size(), "Number of indices to be removed is greater than the size of the vector");

  if (bI.size() > 0)
  {
    Vector u = Vector(v.size() - (int)bI.size());
    int indicesRemoved = 0;
    for (int i = 0; i < v.size(); ++i)
    {
      // Debug
      ASSERT(bI[indicesRemoved] < v.size(), "Index to be removed must be less than the size of the vector");

      if (i != bI[indicesRemoved])
        u[i - indicesRemoved - (i > bI[bI.size() - 1])] = v[i];
      else if (indicesRemoved < bI.size() - 1)
        ++indicesRemoved;
    }
    v = std::move(u);
  }
}

void EquationSystem2D::removeBoundaryIndices(Matrix& A, const std::vector<int> bI) const
{
  // Debug
  ASSERT(A.isSquare(), "Matrix is not square, use other removeBondaryIndices function");
  return removeBoundaryIndices(A, bI, bI);
}

void EquationSystem2D::removeBoundaryIndices(Matrix& A, const std::vector<int> rI, const std::vector<int> cI) const
{
  // Debug
  ASSERT(A.rows() >= rI.size(), "Number of row indices to be removed is greater than the number of rows");
  ASSERT(A.columns() >= cI.size(), "Number of column indices to be removed is greater than the number of columns");

  if (rI.size() > 0 && cI.size() > 0)
  {
    Matrix M = Matrix(A.rows() - (int)rI.size(), A.columns() - (int)cI.size());
    int rowsRemoved = 0;
    for (int i = 0; i < A.rows(); ++i)
    {
      // Debug
      ASSERT(rI[rowsRemoved] < A.rows(), "Row index to be removed must be less than the number of rows");

      if (i != rI[rowsRemoved])
      {
        int columnsRemoved = 0;
        for (int j = 0; j < A.columns(); ++j)
        {
          // Debug
          ASSERT(cI[columnsRemoved] < A.columns(), "Column index to be removed must be less than the number of columns");

          if (j != cI[columnsRemoved])
            M[i - rowsRemoved - (i > rI[rI.size() - 1])][j - columnsRemoved - (j > cI[cI.size() - 1])] = A[i][j];
          else if (columnsRemoved < cI.size() - 1)
            ++columnsRemoved;
        }
      }
      else if (rowsRemoved < rI.size() - 1)
        ++rowsRemoved;
    }
    A = std::move(M);
  }
  else if (rI.size() > 0 && cI.size() == 0)
  {
    Matrix M = Matrix(A.rows() - (int)rI.size(), A.columns());
    int rowsRemoved = 0;
    for (int i = 0; i < A.rows(); ++i)
    {
      // Debug
      ASSERT(rI[rowsRemoved] < A.rows(), "Row index to be removed must be less than the number of rows");

      if (i != rI[rowsRemoved])
        for (int j = 0; j < A.columns(); ++j)
          M[i - rowsRemoved - (i > rI[rI.size() - 1])][j] = A[i][j];
      else if (rowsRemoved < rI.size() - 1)
        ++rowsRemoved;
    }
    A = std::move(M);
  }
  else if (rI.size() == 0 && cI.size() > 0)
  {
    Matrix M = Matrix(A.rows(), A.columns() - (int)cI.size());
    int columnsRemoved = 0;
    for (int j = 0; j < A.columns(); ++j)
    {
      // Debug
      ASSERT(cI[columnsRemoved] < A.columns(), "Column index to be removed must be less than the number of columns");

      if (j != cI[columnsRemoved])
        for (int i = 0; i < A.rows(); ++i)
          M[i][j - columnsRemoved - (j > cI[cI.size() - 1])] = A[i][j];
      else if (columnsRemoved < cI.size() - 1)
        ++columnsRemoved;
    }
    A = std::move(M);
  }
}
