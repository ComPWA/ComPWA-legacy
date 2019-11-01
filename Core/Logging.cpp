// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <array>
#include <chrono>
#include <ctime>

#include "Core/Logging.hpp"

INITIALIZE_EASYLOGGINGPP

namespace ComPWA {

Logging::Logging(std::string level, std::string filename) {
  el::Configurations DefaultConfig;
  // initialize with default values
  DefaultConfig.setToDefault();
  DefaultConfig.setGlobally(el::ConfigurationType::Format,
                            "%datetime [%level] %msg");
  if (filename.empty()) {
    DefaultConfig.setGlobally(el::ConfigurationType::ToFile, "0");
    LOG(INFO) << "Logging to file disabled!";
  } else {
    DefaultConfig.setGlobally(el::ConfigurationType::Filename, filename);
    DefaultConfig.setGlobally(el::ConfigurationType::ToFile, "1");
    LOG(INFO) << "Log file: " << filename;
  }
  // reconfigure the default logger
  el::Loggers::reconfigureLogger("default", DefaultConfig);

  setLogLevel(level);

  LOG(INFO) << "Log level: " << level;

  // Print local time and date at the beginning
  auto time =
      std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

  LOG(INFO) << "Current date and time: " << std::ctime(&time);
};

/// Enable or disable levels TRACE, DEBUG, INFO, WARNING, ERROR, FATAL. An array
/// of strings is passed. E.g {"0","1","1","1","1","1"} to disable TRACE and
/// enable all other levels.
void enableDisableLvl(el::Logger *logger, std::array<std::string, 5> levels) {
  logger->configurations()->set(el::Level::Trace,
                                el::ConfigurationType::Enabled, levels.at(0));
  logger->configurations()->set(el::Level::Debug,
                                el::ConfigurationType::Enabled, levels.at(1));
  logger->configurations()->set(el::Level::Info, el::ConfigurationType::Enabled,
                                levels.at(2));
  logger->configurations()->set(el::Level::Warning,
                                el::ConfigurationType::Enabled, levels.at(3));
  logger->configurations()->set(el::Level::Error,
                                el::ConfigurationType::Enabled, levels.at(4));
}

void Logging::setLogLevel(std::string level) {
  // Capitalize string
  std::transform(level.begin(), level.end(), level.begin(), ::toupper);

  Level = level;

  el::Logger *logger =
      ELPP->registeredLoggers()->get(el::base::consts::kDefaultLoggerId);

  // Normally use the hierarchy mode of easyloggingcpp, e.g.
  // el::Loggers::addFlag(el::LoggingFlag::HierarchicalLogging);
  // el::Loggers::setLoggingLevel(el::Level::Fatal);
  // However, the hierarchy of the easyloggingcpp log levels is currently not
  // convenient and has to be manually reordered to:
  // TRACE, DEBUG, INFO, WARNING, ERROR, FATAL

  if (level == "TRACE")
    enableDisableLvl(logger, {"1", "1", "1", "1", "1"});
  else if (level == "DEBUG")
    enableDisableLvl(logger, {"0", "1", "1", "1", "1"});
  else if (level == "INFO")
    enableDisableLvl(logger, {"0", "0", "1", "1", "1"});
  else if (level == "WARNING")
    enableDisableLvl(logger, {"0", "0", "0", "1", "1"});
  else if (level == "ERROR")
    enableDisableLvl(logger, {"0", "0", "0", "0", "1"});
  else if (level == "FATAL")
    enableDisableLvl(logger, {"0", "0", "0", "0", "0"});
  else {
    enableDisableLvl(logger, {"0", "0", "1", "1", "1"});
    Level = "INFO";
    LOG(WARNING)
        << "Logging::setLogLevel() | Unknown log level " + level +
               ". Available levels are: \"TRACE\", \"DEBUG\", \"INFO\", "
               "\"WARNING\", \"ERROR\", \"FATAL\". Setting log level "
               "to [INFO]!";
  }
  logger->reconfigure();
};

} // namespace ComPWA
