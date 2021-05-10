#pragma once
#include "Precompilied.h"
#include "Utilities/Array2D.h"
#include "Vector.h"

/*
  An nxn matrix of real numbers.
  Indexing starts at 0 and ends at n-1.
*/
class Matrix
{
public:
  Matrix() = delete;

  Matrix(const int size);

  Matrix(const int rows, const int columns);

  Matrix(const Matrix& other) = delete;

  Matrix(Matrix&& other) noexcept;

  Matrix& operator=(const Matrix& other) = delete;

  Matrix& operator=(Matrix&& other) noexcept;

  Container<real>& operator[](const int index);

  const Container<real>& operator[](const int index) const;

  Matrix operator+(const Matrix& other) const;

  void operator*=(const real scalar);

  void operator/=(const real scalar);

  /*
    Removes all entries in the row and column of the specified index
  */
  void removeRowAndCol(const int index);

  /*
    \returns the submatrix resulting from the removal of the first row and specified column
  */
  Matrix subMatrix(const int index) const;

  bool isSquare() const;

  int size() const;

  int rows() const;

  int columns() const;

  void print() const;

  /*
    Very inefficient recursive method to calculate determinant
  */
  real determinant() const;

private:
  int n;
  int m;
  Array2D<real> entries;

  /*
    Recursive determinant helper method
  */
  real determinant(const Matrix& M) const;
};

// --------------------------------------------------------- //

/*
  \returns solution of system Ax = b using LU decomposition.

  Warning: Changes matrix A!
*/
Vector solve(Matrix& A, const Vector& b);

/*
  Passes tempory values to other solve function.

  \returns solution of equation Ax = b.
*/
Vector solve(Matrix&& A, Vector&& b);

Vector solve(std::vector<std::vector<real>>& A, const Vector& b);