#include "Precompilied.h"
#include "Matrix.h"

Matrix::Matrix(const int size)
  : n(size), entries(Array2D<real>(n, n))
{
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      entries[i][j] = 0.0;
}

Matrix::Matrix(Matrix&& other) noexcept
  : n(other.n), entries(std::move(other.entries))
{
}

Container<real>& Matrix::operator[](const int index)
{
  // Debug
  ASSERT(index >= 0, "index must be non-negative");
  ASSERT(index < n, "index must be less than the dimension of the matrix");

  return entries[index];
}

const Container<real>& Matrix::operator[](const int index) const
{
  // Debug
  ASSERT(index >= 0, "index must be non-negative");
  ASSERT(index < n, "index must be less than the dimension of the matrix");

  return entries[index];
}

Matrix Matrix::operator+(const Matrix& other) const
{
  // Debug
  ASSERT(this->size() == other.size(), "Matrices must have the same dimension in order to add");

  Matrix sum = Matrix(n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      sum[i][j] = (*this)[i][j] + other[i][j];
  return sum;
}

void Matrix::operator*=(const real scalar)
{
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      entries[i][j] *= scalar;
}

void Matrix::operator/=(const real scalar)
{
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      entries[i][j] /= scalar;
}

void Matrix::removeRowAndCol(const int index)
{
  Array2D<real> newEntries = Array2D<real>(n - 1, n - 1);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      if (i != index && j != index)
        newEntries[i - (i > index)][j - (j > index)] = entries[i][j];

  entries = std::move(newEntries);
  n--;
}

Matrix Matrix::subMatrix(const int index) const
{
  Matrix subMatrix = Matrix(n - 1);
  for (int i = 1; i < n; ++i)
    for (int j = 0; j < n; ++j)
      if (j != index)
        subMatrix[i - 1][j - (j > index)] = entries[i][j];

  return subMatrix;
}

int Matrix::size() const
{
  return n;
}

void Matrix::print() const
{
  std::cout.precision(16);
  for (int i = 0; i < n; ++i)
  {
    std::cout << "|";
    for (int j = 0; j < n; ++j)
    {
      std::cout << entries[i][j] << " ";
    }
    std::cout << "|" << std::endl;
  }
  std::cout << std::endl;
}

real Matrix::determinant() const
{
  return determinant(*this);
}

real Matrix::determinant(const Matrix& M) const
{
  if (M.size() == 2)
    return M[0][0] * M[1][1] - M[0][1] * M[1][0];
  else
  {
    real sum = 0;
    for (int i = 0; i < M.size(); ++i)
      sum += (i % 2 ? -1 : 1) * M[0][i] * determinant(M.subMatrix(i));
    return sum;
  }
}

Vector solve(Matrix& A, const Vector& b)
{
  ASSERT(A.size() == b.size(), "A and b must be the same dimension");

  const int n = A.size();
  Vector x = Vector(n);

  // Diagonalization
  for (int k = 0; k < n; ++k)
    for (int i = k + 1; i < n; ++i)
    {
      A[i][k] /= A[k][k];
      for (int j = k + 1; j < n; ++j)
        A[i][j] -= A[i][k] * A[k][j];
    }

  // Forward substitution to solve Ld = b
  x[0] = b[0];
  for (int i = 1; i < n; ++i)
  {
    real s = 0;
    for (int j = 0; j < i; ++j)
      s += A[i][j] * x[j];
    x[i] = b[i] - s;
  }

  // Backwards substitution to solve Ux = d
  x[n-1] /= A[n-1][n-1];
  for (int i = n - 2; i >= 0; --i)
  {
    real s = 0;
    for (int j = i + 1; j < n; ++j)
      s += A[i][j] * x[j];
    x[i] = (x[i] - s) / A[i][i];
  }

  return x;
}

Vector solve(Matrix&& A, Vector&& b)
{
  Matrix M = std::move(A);
  Vector v = std::move(b);
  return solve(M, v);
}
