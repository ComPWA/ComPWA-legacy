// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <iomanip>
#include <ctime>

#include "Core/Logging.hpp"

INITIALIZE_EASYLOGGINGPP

namespace ComPWA {

Logging::Logging(std::string lvl) {

  // default logger uses default configurations
  el::Loggers::addFlag(el::LoggingFlag::HierarchicalLogging);

  setLogLevel(lvl);

  LOG(INFO) << "Logging to file disabled!";
  LOG(INFO) << "Log level: " << lvl;

  // Print local time and date at the beginning
  char foo[48];
  std::time_t now = std::time(nullptr); // now stores the current time
  if (0 < strftime(foo, sizeof(foo), "[%c %Z] ", std::localtime(&now)))
    LOG(INFO) << "Current date and time: " << foo;
};

Logging::Logging(std::string out, std::string lvl) {

  // Logging to file
  el::Configurations fileConf;
  fileConf.setToDefault();
  fileConf.setGlobally(el::ConfigurationType::Filename, out);
  fileConf.setGlobally(el::ConfigurationType::Format,
                       "%datetime [%level] %msg");
  fileConf.setGlobally(el::ConfigurationType::ToFile, "1");
  fileConf.setGlobally(el::ConfigurationType::ToStandardOutput, "1");
  // default logger uses default configurations
  el::Loggers::reconfigureLogger("default", fileConf);

  el::Loggers::addFlag(el::LoggingFlag::HierarchicalLogging);

  setLogLevel(lvl);


  LOG(INFO) << "Log file: " << out;
  LOG(INFO) << "Log level: " << lvl;

  // Print local time and date at the beginning
  
  // GCC < 5.0: std::put_time is not available
  // std::time_t now = std::time(nullptr);
  //  LOG(INFO) << "Current date and time: "
  //            << std::put_time(std::localtime(&now), "%c %Z");

  char foo[48];
  std::time_t now = std::time(nullptr); // now stores the current time
  if (0 < strftime(foo, sizeof(foo), "[%c %Z] ", std::localtime(&now)))
    LOG(INFO) << "Current date and time: " << foo;
};

void Logging::setLogLevel(std::string minLevel) {

  // Capitalize string
  std::transform(minLevel.begin(), minLevel.end(), minLevel.begin(), ::toupper);

  if (minLevel == "TRACE")
    el::Loggers::setLoggingLevel(el::Level::Trace);
  else if (minLevel == "DEBUG")
    el::Loggers::setLoggingLevel(el::Level::Debug);
  else if (minLevel == "FATAL")
    el::Loggers::setLoggingLevel(el::Level::Fatal);
  else if (minLevel == "ERROR")
    el::Loggers::setLoggingLevel(el::Level::Error);
  else if (minLevel == "WARNING")
    el::Loggers::setLoggingLevel(el::Level::Warning);
  else if (minLevel == "INFO")
    el::Loggers::setLoggingLevel(el::Level::Info);
  else
    throw std::runtime_error("Logging::setLogLevel() | Log level " + minLevel +
                             " unknown.");
};

} // namespace ComPWA
