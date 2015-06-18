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

#include "HelicityDecayXMLConfigReader.hpp"

#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>

namespace HelicityFormalism {

using boost::property_tree::ptree;

HelicityDecayXMLConfigReader::HelicityDecayXMLConfigReader(
    HelicityDecayConfiguration &decay_configuration) :
    decay_configuration_(decay_configuration) {
}

HelicityDecayXMLConfigReader::~HelicityDecayXMLConfigReader() {
}

std::vector<unsigned int> HelicityDecayXMLConfigReader::parseIDList(
    const boost::property_tree::ptree &pt) const {
  std::vector<unsigned int> ids;
  BOOST_FOREACH(ptree::value_type const &id_list_item, pt.get_child("id_list")){
  ids.push_back(id_list_item.second.get<unsigned int>(""));
}
  return ids;
}

void HelicityDecayXMLConfigReader::readConfig(const std::string &filename) {
  // Create an empty property tree object

  // Load the XML file into the property tree. If reading fails
  // (cannot open file, parse error), an exception is thrown.
  read_xml(filename, pt_, boost::property_tree::xml_parser::trim_whitespace);
  // Get the filename and store it in the m_file variable.
  // Note that we construct the path to the value by separating
  // the individual keys with dots. If dots appear in the keys,
  // a path type with a different separator can be used.
  // If the amplitude_setup.filename key is not found, an
  // exception is thrown.
  pt_ = pt_.get_child("amplitude_setup");

  // Iterate over the amplitude_setup.resonances section and
  // store all found resonance names in the m_name set.
  // The get_child() function returns a reference to the child
  // at the specified path; if there is no such child, it throws.
  // Property tree iterators are models of BidirectionalIterator.
  BOOST_FOREACH(ptree::value_type const& v, pt_.get_child("FinalState")){
  ParticleState ps;
  ps.J_ = 0;
  ps.M_ = 0;
  ps.particle_id_ = v.second.get<unsigned int>("id");
  ps.name_ = v.second.get<std::string>("name");
  decay_configuration_.addFinalStateParticle(ps);
}
  BOOST_FOREACH(ptree::value_type const& v, pt_.get_child("IntermediateStates")){
  ParticleState ps;
  ps.J_ = 0;
  ps.M_ = 0;
  ps.particle_id_ = v.second.get<unsigned int>("id");
  ps.name_ = v.second.get<std::string>("name");
  decay_configuration_.addFinalStateParticle(ps);
  decay_configuration_.addIntermediateStateParticle(ps);
}

  BOOST_FOREACH(ptree::value_type const& decay_tree, pt_.get_child("DecayTrees")){
  BOOST_FOREACH(ptree::value_type const& two_body_decay, decay_tree.second) {
    std::vector<unsigned int> mother_ids = parseIDList(two_body_decay.second.get_child("particle_ids"));
    std::vector<unsigned int> daughter_1_ids = parseIDList(two_body_decay.second.get_child("daughter_1_ids"));
    std::vector<unsigned int> daughter_2_ids = parseIDList(two_body_decay.second.get_child("daughter_2_ids"));
    std::pair<std::vector<unsigned int>, std::vector<unsigned int> > daughter_ids = std::make_pair(daughter_1_ids, daughter_2_ids);
    for(unsigned int i = 0; i < mother_ids.size(); ++i) {
      decay_configuration_.addDecayToCurrentDecayTopology(mother_ids[i], daughter_ids);
    }
  }
  decay_configuration_.addCurrentDecayTopologyToList();
}

}

}
/* namespace HelicityFormalism */
