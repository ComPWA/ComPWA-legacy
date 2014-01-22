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
	Logging(std::string outFileName="output.log"){ init(outFileName); };
private:
	void init(std::string out){
		boost::log::add_common_attributes();
		boost::log::add_file_log( boost::log::keywords::file_name=out,
				//			keywords::format="(%LineID%) [%TimeStamp%][%Severity%]: %Message%"
				boost::log::keywords::format =
						(
								boost::log::expressions::stream
								<< boost::log::expressions::format_date_time< boost::posix_time::ptime >("TimeStamp", "%Y-%m-%d %H:%M:%S")
								<< " [" << boost::log::trivial::severity<< "] : "
								<< boost::log::expressions::smessage
						)
		);
		boost::log::add_console_log(std::cout,
				boost::log::keywords::format =
						(
								boost::log::expressions::stream
								<< boost::log::expressions::format_date_time< boost::posix_time::ptime >("TimeStamp", "%H:%M:%S")
								<< " [" << std::setw(7) << boost::log::trivial::severity<< "] : "
								<< boost::log::expressions::smessage
						)
		);
		boost::log::core::get()->set_filter(boost::log::trivial::severity >= boost::log::trivial::debug);
	}
};



#endif /* LOGGING_HPP_ */
