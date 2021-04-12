#pragma once
#include "StdLib.h"
#include "Log.h"

using real = double;
using real1DFunction = real(*)(real);
using real2DFunction = real(*)(real, real);

// Important constants
constexpr real PI = 3.14159265358979323846;
constexpr real TOLERANCE = 1e-14;


// ------------------------------ Macros ------------------------------ //

/*
  Gets file name of the file in which this macro is invoked.
  Found here: https://stackoverflow.com/questions/8487986/file-macro-shows-full-path
*/
#define __FILENAME__ (strrchr(__FILE__, '\\') ? strrchr(__FILE__, '\\') + 1 : __FILE__)

/*
  Assert macro that ends program if condition is not met.
  Will print out filename and line where the error occured.

  The assert statments will ONLY appear in debug mode, and
  will be replaced with empty lines in release mode.
*/
#if FE_DEBUG
#define ASSERT(condition, messageOtherwise) if (!(condition)) logError(messageOtherwise, LogLevel::Error, __FILENAME__, __LINE__)
#else
#define ASSERT(condition, messageOtherwise)
#endif

/*
  Logs a message to the console, along with the file name and line
  where the LOG macro is invoked.  Uses the LogLevel enum class as a second argument.
  There are three log levels:

  Error - Meant for when things go seriously wrong.  Terminates program immediately.
  Warning - Meant for when something worrying or unexpected has happened, but is not fatal to the program.
  Info - Meant to merely provide additional information.
*/
#define LOG(message, logLevel) logError(message, logLevel, __FILENAME__, __LINE__)