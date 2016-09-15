/*
 * DecayGeneratorConfig.cpp
 *
 *  Created on: Aug 22, 2016
 *      Author: steve
 */

#include <boost/property_tree/json_parser.hpp>
#include <boost/log/trivial.hpp>

#include "Physics/DecayTree/DecayGeneratorConfig.hpp"

using boost::property_tree::ptree;

namespace ComPWA {
namespace Physics {
namespace DecayTree {

DecayGeneratorConfig::DecayGeneratorConfig() {
}

DecayGeneratorConfig::~DecayGeneratorConfig() {
}

void DecayGeneratorConfig::readConfig(const std::string& config_url) {
  BOOST_LOG_TRIVIAL(info)
      << "DecayGeneratorConfig::readConfig(): reading decay generator config file "
      << config_url;
  read_json(config_url, config_ptree_);
}

const boost::property_tree::ptree& DecayGeneratorConfig::getConfig() const {
  return config_ptree_;
}

} /* namespace DecayTree */
} /* namespace Physics */
} /* namespace ComPWA */
