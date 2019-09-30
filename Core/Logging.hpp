// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef LOGGING_HPP_
#define LOGGING_HPP_

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
  Logging(std::string level = "INFO", std::string filename = "");

  void setLogLevel(std::string level);

  std::string getLogLevel() { return Level; };

private:
  std::string Level;
};

} // namespace ComPWA

#endif
