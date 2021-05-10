#include "Precompilied.h"
#include "Vector.h"

Vector::Vector(const int size)
  : n(size)
{
  entries = new real[n];
  for (int i = 0; i < n; ++i)
    entries[i] = 0.0;
}

Vector::Vector(Vector&& other) noexcept
  : n(other.n)
{
  entries = other.entries;
  other.entries = nullptr;
}

Vector& Vector::operator=(Vector&& other) noexcept
{
  if (&other != this)
  {
    n = other.n;
    delete[] entries;
    entries = other.entries;
    other.entries = nullptr;
  }
  return *this;
}

real& Vector::operator[](const int index)
{
  // Debug
  ASSERT(index >= 0, "index must be non-negative");
  ASSERT(index < n, "index must be less than the dimension of the vector");

  return entries[index];
}

const real& Vector::operator[](const int index) const
{
  // Debug
  ASSERT(index >= 0, "index must be non-negative");
  ASSERT(index < n, "index must be less than the dimension of the vector");

  return entries[index];
}

Vector Vector::operator+(const Vector& other) const
{
  // Debug
  ASSERT(this->size() == other.size(), "Vectors must have the same number of elements in order to add");

  Vector sum = Vector(n);
  for (int i = 0; i < n; ++i)
    sum[i] = (*this)[i] + other[i];
  return sum;
}

void Vector::remove(const int index)
{
  // Debug
  ASSERT(index >= 0, "index must be non-negative");
  ASSERT(index < n, "index must be less than the dimension of the vector");

  real* newEntries = new real[n - 1];
  for (int i = 0; i < n; ++i)
    if (i != index)
      newEntries[i - (i > index)] = entries[i];

  delete[] entries;
  entries = newEntries;
  n--;
}

void Vector::insert(const real a, const int index)
{
  // Debug
  ASSERT(index >= 0, "index must be non-negative");
  ASSERT(index <= n, "index must be less than or equal to the dimension of the vector");

  real* newEntries = new real[n + 1];
  for (int i = 0; i < n; ++i)
  {
    if (i == index)
      newEntries[i] = a;
    newEntries[i + (i >= index)] = entries[i];
  }
  if (index == n)
    newEntries[n] = a;

  delete[] entries;
  entries = newEntries;
  n++;
}

void Vector::print() const
{
  std::cout.precision(16);
  for (int i = 0; i < n; ++i)
  {
    std::cout << "|" << entries[i] << " " << "|" << std::endl;
  }
  std::cout << std::endl;
}

int Vector::size() const
{
  return n;
}
