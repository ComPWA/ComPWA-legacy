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

#ifndef HELICITYDECAYXMLCONFIGREADER_HPP_
#define HELICITYDECAYXMLCONFIGREADER_HPP_

#include "HelicityDecayConfiguration.hpp"

#include <boost/property_tree/ptree.hpp>

namespace HelicityFormalism {

class HelicityDecayXMLConfigReader {
  HelicityDecayConfiguration &decay_configuration_;

  boost::property_tree::ptree pt_;

  std::vector<unsigned int> parseIDList(
      const boost::property_tree::ptree &pt) const;

public:
  HelicityDecayXMLConfigReader(HelicityDecayConfiguration &decay_configuration);
  virtual ~HelicityDecayXMLConfigReader();

  void readConfig(const std::string &filename);
};

} /* namespace HelicityFormalism */

#endif /* HELICITYDECAYXMLCONFIGREADER_HPP_ */
