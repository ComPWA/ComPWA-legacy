/*
 * DecayGeneratorConfig.hpp
 *
 *  Created on: Aug 22, 2016
 *      Author: steve
 */

#ifndef PHYSICS_DECAYTREE_DECAYGENERATORCONFIG_HPP_
#define PHYSICS_DECAYTREE_DECAYGENERATORCONFIG_HPP_

#include <boost/property_tree/ptree.hpp>

namespace ComPWA {
namespace Physics {
namespace DecayTree {

class DecayGeneratorConfig {
  boost::property_tree::ptree config_ptree_;

  DecayGeneratorConfig();
public:
  static DecayGeneratorConfig& Instance() {
    static DecayGeneratorConfig instance;
    return instance;
  }
  virtual ~DecayGeneratorConfig();

  DecayGeneratorConfig(DecayGeneratorConfig const&) = delete;
  void operator=(DecayGeneratorConfig const&) = delete;

  void readConfig(const std::string& config_url);

  const boost::property_tree::ptree& getConfig() const;
};

} /* namespace DecayTree */
} /* namespace Physics */
} /* namespace ComPWA */

#endif /* PHYSICS_DECAYTREE_DECAYGENERATORCONFIG_HPP_ */
