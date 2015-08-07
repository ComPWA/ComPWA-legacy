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

#include "DecayXMLConfigReader.hpp"

#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>

namespace HelicityFormalism {

using boost::property_tree::ptree;

DecayXMLConfigReader::DecayXMLConfigReader(
    DecayConfiguration &decay_configuration) :
    decay_configuration_(decay_configuration) {
}

DecayXMLConfigReader::~DecayXMLConfigReader() {
}

std::vector<unsigned int> DecayXMLConfigReader::parseIDList(
    const boost::property_tree::ptree &pt) const {
  std::vector<unsigned int> ids;
  BOOST_FOREACH(ptree::value_type const &id_list_item, pt.get_child("id_list")){
    ids.push_back(id_list_item.second.get<unsigned int>(""));
  }
  return ids;
}

void DecayXMLConfigReader::readConfig(const std::string &filename) {
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
  pt_ = pt_.get_child("DecayTreeSetup");

  // Iterate over the amplitude_setup.resonances section and
  // store all found resonance names in the m_name set.
  // The get_child() function returns a reference to the child
  // at the specified path; if there is no such child, it throws.
  // Property tree iterators are models of BidirectionalIterator.
  BOOST_FOREACH(ptree::value_type const& v, pt_.get_child("FinalState")) {
    ParticleState ps;
    ps.helicity_state_information.J_ = 0;
    ps.helicity_state_information.M_ = 0;
    ps.id_ = v.second.get<unsigned int>("id");
    ps.name_ = v.second.get<std::string>("name");
    decay_configuration_.addFinalStateParticle(ps);
  }
  BOOST_FOREACH(ptree::value_type const& v, pt_.get_child("IntermediateStates")) {
    ParticleState ps;
    ps.helicity_state_information.J_ = 0;
    ps.helicity_state_information.M_ = 0;
    ps.id_ = v.second.get<unsigned int>("id");
    ps.name_ = v.second.get<std::string>("name");
    decay_configuration_.addFinalStateParticle(ps);
    decay_configuration_.addIntermediateStateParticle(ps);
  }

  BOOST_FOREACH(ptree::value_type const& decay_tree, pt_.get_child("DecayTrees")) {
    BOOST_FOREACH(ptree::value_type const& decay_node, decay_tree.second) {
      std::vector<unsigned int> mother_ids = parseIDList(decay_node.second.get_child("mother"));
      std::vector<std::vector<unsigned int> > daughters_id_lists;
      BOOST_FOREACH(ptree::value_type const& daugthers, decay_node.second.get_child("daughters")) {
        daughters_id_lists.push_back(parseIDList(daugthers.second.get_child("daughter")));
      }
      for(unsigned int i = 0; i < mother_ids.size(); ++i) {
        decay_configuration_.addDecayToCurrentDecayTree(mother_ids[i], daughters_id_lists);
      }
    }
    decay_configuration_.addCurrentDecayTreeToList();
  }
}

}
/* namespace HelicityFormalism */
