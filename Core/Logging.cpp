/*
 * Logging.cpp
 *
 *  Created on: Mar 19, 2014
 *      Author: weidenka
 */


#include "Core/Logging.hpp"

void Logging::init(std::string out,boost::log::trivial::severity_level minLevel){
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
	boost::log::core::get()->set_filter(boost::log::trivial::severity >= minLevel);
	BOOST_LOG_TRIVIAL(info)<<"Logging: using output filename: "<<out<<", Severity level: "<<minLevel;
}
