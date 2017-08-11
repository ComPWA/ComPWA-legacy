// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Core/Logging.hpp"

#include <boost/log/expressions.hpp>
#include <boost/log/expressions/formatters/date_time.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/log/support/date_time.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/date_time/posix_time/ptime.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

namespace ComPWA {
using namespace boost::log;
  
void Logging::init(std::string out, trivial::severity_level minLevel) {
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
    LOG(info)
        << "Logging::init() | Logging to file disabled. Console severity level: "
        << minLevel;
  } else {
    add_file_log(
        keywords::file_name = out,
        //			keywords::format="(%LineID%)
        //[%TimeStamp%][%Severity%]: %Message%"
        keywords::format =
            (expressions::stream
             << expressions::format_date_time<boost::posix_time::ptime>(
                    "TimeStamp", "%Y-%m-%d %H:%M:%S")
             << " [" << trivial::severity << "] : " << expressions::smessage));
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
}

void Logging::SetLogLevel(trivial::severity_level minLevel) {
  core::get()->set_filter(trivial::severity >= minLevel);
  LOG(info) << "New severity level: " << minLevel;
}

} /* namespace ComPWA */
