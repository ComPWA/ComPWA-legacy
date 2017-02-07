/*
 * Logging.hpp
 *
 *  Created on: Jan 21, 2014
 *      Author: weidenka
 */

#ifndef LOGGING_HPP_
#define LOGGING_HPP_


#include <boost/log/trivial.hpp>
#include <boost/log/common.hpp>


namespace ComPWA {
    
    //Redefine BOOST_LOG_TRIVIAL(lvl) to LOG(lvl
#define LOG(lvl) \
    BOOST_LOG_STREAM_WITH_PARAMS(::boost::log::trivial::logger::get(),\
(::boost::log::keywords::severity = ::boost::log::trivial::lvl))
    
    class Logging
    {
    public:
        Logging(std::string outFileName="output.log",
                boost::log::trivial::severity_level minLevel
                = boost::log::trivial::debug){
            init(outFileName,minLevel);
        };
        
        void setLogLevel(boost::log::trivial::severity_level minLevel);
        
    private:
        void init(std::string out,boost::log::trivial::severity_level minLevel);
    };
    
} /* namespace ComPWA */

#endif /* LOGGING_HPP_ */
