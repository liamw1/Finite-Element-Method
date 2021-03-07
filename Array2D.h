#pragma once
#include "Precompilied.h"
#include "Container.h"

/*
  A 2D-style array that stores its members in a single
  heap-allocated block in memory.

  Elements can be accessed with brackets:
  arr[i][j]

  Huge advantages in efficiency for "tall" arrays (rows >> cols).
*/
template<typename T>
class Array2D
{
public:
  Array2D() = delete;

  Array2D(const int size1, const int size2)
    : container(Container<T>(size1, size2))
  {
  }

  Array2D(const Array2D& other) = delete;

  Array2D(Array2D&& other) noexcept
    : container(std::move(other.container))
  {
  }

  Array2D& operator=(const Array2D& other) = delete;

  Array2D& operator=(Array2D&& other) noexcept
  {
    if (&other != this)
      container = std::move(other.container);
    return *this;
  }

  Container<T>& operator[](const int index)
  {
    container.index1 = index;
    return container;
  }

  const Container<T>& operator[](const int index) const
  {
    container.index1 = index;
    return container;
  }

private:
  Container<T> container;
};
