// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <chrono>
#include <boost/log/support/date_time.hpp>
#include "Core/Logging.hpp"

INITIALIZE_EASYLOGGINGPP

namespace ComPWA {
using namespace boost::log;

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

  // Logging to terminal
  //  el::Logger* terminalLog = el::Loggers::getLogger("terminalLog");
  //  el::Configurations terminalConf;
  //  terminalConf.setToDefault();
  //  terminalConf.setGlobally(el::ConfigurationType::Format,
  //                       "%datetime [%level] %msg");
  //  terminalConf.setGlobally(el::ConfigurationType::ToFile, "0");
  //  terminalConf.setGlobally(el::ConfigurationType::ToStandardOutput, "1");
  //  el::Loggers::reconfigureLogger("terminalLog", terminalConf);

  // Print local time and date at the beginning
  boost::posix_time::ptime todayUtc(
      boost::gregorian::day_clock::universal_day(),
      boost::posix_time::second_clock::local_time().time_of_day());
  LOG(INFO) << "Log file: " << out;
  LOG(INFO) << "Current date and time: "
            << boost::posix_time::to_simple_string(todayUtc);
};

void Logging::setLogLevel(std::string minLevel){
    // ToDo: implement setLogLevel
};

} // namespace ComPWA
