// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef LOGGING_HPP_
#define LOGGING_HPP_

// #include "easylogging++.h"
#include "ThirdParty/easyloggingpp/easylogging++.h"

namespace ComPWA {

///
/// \class Logging
/// Logging class provides an interface for logging all over the framework.
/// Behind the scenes easyloggingcpp is currently used which allows a detailed
/// on logging format and log levels
///

class Logging {
public:
  /// Logging to file and stdout with level minLevel
  Logging(std::string outFileName,
          std::string minLevel = "DEBUG");

  /// Logging to stdout only with level minLevel
  Logging(std::string minLevel = "DEBUG");

  void setLogLevel(std::string minLevel);
};

} // namespace ComPWA

#endif
