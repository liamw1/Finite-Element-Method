#include "Precompilied.h"
#include "HomeworkDrivers.h"

std::ostream& operator<<(std::ostream& os, const BC_Type BC)
{
  os << (int)BC;
  return os;
}

int main()
{
  Hwk10_C1_Driver();
}