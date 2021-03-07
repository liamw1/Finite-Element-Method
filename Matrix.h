#pragma once
#include "Precompilied.h"
#include "Array2D.h"
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

  Matrix(const Matrix& other) = delete;

  Matrix(Matrix&& other) noexcept;

  Matrix& operator=(const Matrix& other) = delete;

  Container<real>& operator[](const int index);

  const Container<real>& operator[](const int index) const;

  Matrix operator+(const Matrix& other) const;

  /*
    Removes all entries in the row and column of the specified index
  */
  void removeRowAndCol(const int index);

  int size() const;

  void print() const;

private:
  int n;
  Array2D<real> entries;
};

// --------------------------------------------------------- //

/*
  Solves the linear system Ax = b using LU decomposition.
  Returns the vector x.

  Warning: Changes matrix A!
*/
Vector solve(Matrix& A, const Vector& b);

/*
  Solves the linear system Ax = b using LU decomposition.
  Returns the vector x.
*/
Vector solve(Matrix&& A, Vector&& b);