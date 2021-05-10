#pragma once
#include "Precompilied.h"

struct Timer
{
  std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
  std::chrono::duration<float> duration;
  const char* timerName;

  Timer() = delete;

  Timer(const char* name);

  Timer(const Timer& other) = delete;

  Timer(Timer&& other) noexcept;

  Timer& operator=(Timer& other) = delete;

  void timeStart();

  void timeEnd();

  ~Timer();
};