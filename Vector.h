#pragma once
#include "Precompilied.h"

/*
  An n-dimensional vector of real numbers.
  Indexing starts at 0 and ends at n-1.
*/
class Vector
{
public:
  Vector() = delete;

  Vector(const int size);

  Vector(const Vector& other) = delete;

  Vector(Vector&& other) noexcept;

  Vector& operator=(Vector& other) = delete;

  Vector& operator=(Vector&& other) noexcept;

  /*
    Accesses element of the vector at the specified index.
  */
  real& operator[](const int index);

  /*
    Accesses element of the vector at the specified index.
  */
  const real& operator[](const int index) const;

  Vector operator+(const Vector& other) const;

  /*
    Removes element at specified index.  
    Vector will be one less in size.
  */
  void remove(const int index);

  /*
    Inserts a number "a" at the specified index.  
    Vector will be one greater in size.
  */
  void insert(const real a, const int index);

  /*
    Prints the vector to the console as a column vector.
  */
  void print() const;

  /*
    Returns size of the vector.
  */
  int size() const;

private:
  real* entries = nullptr;
  int n;
};