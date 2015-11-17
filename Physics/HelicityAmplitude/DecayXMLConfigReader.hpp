//-------------------------------------------------------------------------------
// Copyright (c) 2013 Stefan Pflueger.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//   Stefan Pflueger - initial API and implementation
//-------------------------------------------------------------------------------

#ifndef PHYSICS_HELICITYAMPLITUDE_DECAYXMLCONFIGREADER_HPP_
#define PHYSICS_HELICITYAMPLITUDE_DECAYXMLCONFIGREADER_HPP_

#include <boost/property_tree/ptree.hpp>

#include "DecayConfiguration.hpp"

namespace HelicityFormalism {

class DecayXMLConfigReader {
  DecayConfiguration &decay_configuration_;

  boost::property_tree::ptree pt_;

  std::map<unsigned int, ParticleStateInfo> template_particle_states_;

  ParticleStateInfo parseParticleStateBasics(
      const boost::property_tree::ptree &pt) const;

  ParticleStateInfo parseParticleStateRemainders(
      const boost::property_tree::ptree &pt) const;

  //std::vector<unsigned int> parseIDList(
  //    const boost::property_tree::ptree &pt) const;

public:
  DecayXMLConfigReader(DecayConfiguration &decay_configuration);
  virtual ~DecayXMLConfigReader();

  void readConfig(const std::string &filename);
};

} /* namespace HelicityFormalism */

#endif /* PHYSICS_HELICITYAMPLITUDE_DECAYXMLCONFIGREADER_HPP_ */
