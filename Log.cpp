#include "Log.h"

static std::vector<std::pair<std::string, int>> errors;

void logError(std::string message, LogLevel logLevel, std::string file, int line)
{
  bool uniqueError = true;
  for (int i = 0; i < errors.size(); ++i)
  {
    if (line == errors[i].second && file == errors[i].first)
    {
      uniqueError = false;
      break;
    }
  }

  if (uniqueError)
  {
    switch (logLevel)
    {
    case LogLevel::Error:
      std::cout << "Error: ";
      break;
    case LogLevel::Warning:
      std::cout << "Warning: ";
      break;
    case LogLevel::Info:
      std::cout << "Info: ";
      break;
    }

    std::pair<std::string, int> error = { file, line };
    errors.emplace_back(error);
    std::cout << message << " ---> " << file << " line " << line << std::endl;

    if (logLevel == LogLevel::Error)
      exit(1);
  }
}