#pragma once
#include "StdLib.h"

enum class LogLevel
{
  Error,
  Warning,
  Info
};

void logError(std::string message, LogLevel logLevel, std::string file, int line);