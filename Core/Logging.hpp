/*
 * Logging.hpp
 *
 *  Created on: Jan 21, 2014
 *      Author: weidenka
 */

#ifndef LOGGING_HPP_
#define LOGGING_HPP_

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include "boost/log/expressions/formatters/date_time.hpp"
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/log/support/date_time.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/date_time/posix_time/ptime.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

//using namespace boost::log;
class Logging
{
public:
	Logging(std::string outFileName="output.log",
			boost::log::trivial::severity_level minLevel=boost::log::trivial::debug){
		init(outFileName,minLevel);
	};
	void setLogLevel(boost::log::trivial::severity_level minLevel);
private:
	void init(std::string out,boost::log::trivial::severity_level minLevel);
};
#endif /* LOGGING_HPP_ */
