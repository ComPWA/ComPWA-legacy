// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef LOGGING_HPP_
#define LOGGING_HPP_

#include <boost/log/trivial.hpp>
#include <boost/log/common.hpp>

namespace ComPWA {

// Redefine LOG(lvl) to LOG(lvl
#define LOG(lvl)                                                               \
  BOOST_LOG_STREAM_WITH_PARAMS(                                                \
      ::boost::log::trivial::logger::get(),                                    \
      (::boost::log::keywords::severity = ::boost::log::trivial::lvl))

inline boost::log::trivial::severity_level
stringToLoggingLevel(std::string logStr) {
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

///
/// \class Logging
/// Logging class privides an interface for logging all over the framework.
/// Behind the scenes boost::log is currently used which allows a detailed
/// on logging format and log levels
///
class Logging {
public:
  Logging(std::string outFileName = "output.log",
          boost::log::trivial::severity_level minLevel =
              boost::log::trivial::debug) {
    init(outFileName, minLevel);
  };

  void setLogLevel(boost::log::trivial::severity_level minLevel);

  void setLogLevel(std::string minLevel) {
    setLogLevel(stringToLoggingLevel(minLevel));
  };

private:
  void init(std::string out, boost::log::trivial::severity_level minLevel);
};

} // namespace ComPWA

#endif
