// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Core/Logging.hpp"

#include <boost/log/trivial.hpp>
#include <boost/log/common.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/log/support/date_time.hpp>

namespace ComPWA {
using namespace boost::log;

boost::log::trivial::severity_level stringToLoggingLevel(std::string logStr) {
  boost::log::trivial::severity_level logLv;
  if (logStr == "trace")
    logLv = boost::log::trivial::trace;
  else if (logStr == "debug")
    logLv = boost::log::trivial::debug;
  else if (logStr == "info")
    logLv = boost::log::trivial::info;
  else if (logStr == "warning")
    logLv = boost::log::trivial::warning;
  else if (logStr == "error")
    logLv = boost::log::trivial::error;
  else if (logStr == "fatal")
    logLv = boost::log::trivial::fatal;
  else
    throw std::runtime_error("Logging stringToLogLevel | unknown log level \"" +
                             logStr + "\"");
  return logLv;
}

Logging::Logging(std::string out, std::string lvl) {
  auto minLevel = stringToLoggingLevel(lvl);
  add_common_attributes();
  boost::log::add_console_log(
      std::cout,
      keywords::format =
          (expressions::stream
           << expressions::format_date_time<boost::posix_time::ptime>(
                  "TimeStamp", "%H:%M:%S")
           << " [" << std::setw(7) << trivial::severity
           << "] : " << expressions::smessage));
  if (out == "") {
    LOG(info) << "Logging::init() | Logging to file disabled. Console severity "
                 "level: "
              << minLevel;
  } else {
    add_file_log(
        keywords::file_name = out,
        //      keywords::format="(%LineID%)
        //[%TimeStamp%][%Severity%]: %Message%"
        keywords::format =
            (expressions::stream
             << expressions::format_date_time<boost::posix_time::ptime>(
                    "TimeStamp", "%Y-%m-%d %H:%M:%S")
             << " [" << std::setw(7) << trivial::severity
             << "] : " << expressions::smessage));
    LOG(info) << "Logging: using output filename: " << out
              << ", Severity level: " << minLevel;
  }
  core::get()->set_filter(trivial::severity >= minLevel);

  // Print local time and date at the beginning
  boost::posix_time::ptime todayUtc(
      boost::gregorian::day_clock::universal_day(),
      boost::posix_time::second_clock::local_time().time_of_day());
  LOG(info) << "Current date and time: "
            << boost::posix_time::to_simple_string(todayUtc);
};

void Logging::log(logLvl level){
  BOOST_LOG_STREAM_WITH_PARAMS(
      ::boost::log::trivial::logger::get(),
      (::boost::log::keywords::severity = ::boost::log::trivial::debug));
}

void Logging::setLogLevel(std::string minLevel) {
  core::get()->set_filter(trivial::severity >= stringToLoggingLevel(minLevel));
  LOG(info) << "New severity level: " << minLevel;
};

} // namespace ComPWA
