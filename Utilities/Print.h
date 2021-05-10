#pragma once
#include "Precompilied.h"

void print();

/*
  A print function with variable number of arguments.
*/
template<typename T, typename... Types>
void print(const T& head, const Types&... tail)
{
  std::cout.precision(16);
  std::cout << head;
  print(tail...);
}