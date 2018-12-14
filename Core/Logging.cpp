// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <chrono>
#include <ctime>

#include "Core/Logging.hpp"

INITIALIZE_EASYLOGGINGPP

namespace ComPWA {

Logging::Logging(std::string out, std::string lvl) {

  // Logging to file
  el::Configurations DefaultConfig;
  DefaultConfig.setToDefault();
  DefaultConfig.setGlobally(el::ConfigurationType::Filename, out);
  DefaultConfig.setGlobally(el::ConfigurationType::Format,
                            "%datetime [%level] %msg");
  DefaultConfig.setGlobally(el::ConfigurationType::ToFile, "1");
  DefaultConfig.setGlobally(el::ConfigurationType::ToStandardOutput, "1");

  // default logger uses default configurations
  el::Loggers::reconfigureLogger("default", DefaultConfig);

  setLogLevel(DefaultConfig, lvl);

  LOG(INFO) << "Log file: " << out;
  LOG(INFO) << "Log level: " << lvl;

  // Print local time and date at the beginning
  auto time =
      std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

  LOG(INFO) << "Current date and time: " << std::ctime(&time);
};

void Logging::setLogLevel(el::Configurations &Config, std::string Level) const {
  // Capitalize string
  std::transform(Level.begin(), Level.end(), Level.begin(), ::toupper);

  // Normally use the hierarchy mode of easyloggingcpp, e.g.
  // el::Loggers::addFlag(el::LoggingFlag::HierarchicalLogging);
  // el::Loggers::setLoggingLevel(el::Level::Fatal);
  // However, the hierarchy of the easyloggingcpp log levels is currently not
  // convenient and has to be manually reordered to:
  // TRACE, DEBUG, INFO, WARNING, ERROR, FATAL

  std::vector<el::Level> OffLevels;
  if (Level == "TRACE") {
    // trace is the highest level, outputs everything
  } else if (Level == "DEBUG")
    OffLevels = {el::Level::Trace};
  else if (Level == "INFO")
    OffLevels = {el::Level::Trace, el::Level::Debug};
  else if (Level == "WARNING")
    OffLevels = {el::Level::Trace, el::Level::Debug, el::Level::Info};
  else if (Level == "ERROR")
    OffLevels = {el::Level::Trace, el::Level::Debug, el::Level::Info,
                 el::Level::Warning};
  else if (Level == "FATAL")
    OffLevels = {el::Level::Trace, el::Level::Debug, el::Level::Info,
                 el::Level::Warning, el::Level::Error};
  else {
    OffLevels = {el::Level::Trace, el::Level::Debug};
    LOG(WARNING) << "Logging::setLogLevel() | Log level " + Level +
                        " unknown. Setting log level to [Info] instead!";
  }
  for (auto x : OffLevels) {
    Config.set(x, el::ConfigurationType::Enabled, "0");
  }
  el::Loggers::reconfigureLogger("default", Config);
};

} // namespace ComPWA
