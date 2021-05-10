#include "Precompilied.h"
#include "Timer.h"

Timer::Timer(const char* name)
  : timerName(name)
{
  duration = std::chrono::high_resolution_clock::duration::zero();
}

Timer::Timer(Timer&& other) noexcept
  : start(std::move(other.start)), end(std::move(other.end)), duration(std::move(other.duration))
{
  timerName = other.timerName;
  other.timerName = nullptr;
}

void Timer::timeStart()
{
  start = std::chrono::high_resolution_clock::now();
}

void Timer::timeEnd()
{
  end = std::chrono::high_resolution_clock::now();
  duration += end - start;
}

Timer::~Timer()
{
  std::cout << timerName << " took " << duration.count() * 1000.0f << "ms" << std::endl;
}