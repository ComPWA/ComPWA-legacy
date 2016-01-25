/*
 * PhysConst.cpp
 *
 *  Created on: Oct 22, 2013
 *      Author: weidenka
 */

#include <stdlib.h>
#include <sstream>
#include <string>
#include <exception>

#include "Core/PhysConst.hpp"
// Boost header files go here
#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/log/trivial.hpp>
using namespace boost::log;

namespace ComPWA {

PhysConst::PhysConst() {
  const char* pPath = getenv("COMPWA_DIR");
  std::string path = "";
  try {
    path = std::string(pPath);
  }
  catch (std::logic_error) {
    BOOST_LOG_TRIVIAL(error)<<"Environment Variable COMPWA_DIR not set?"<<std::endl;
  }
  particleFileName = path + "/Physics/particles.xml";
  particleDefaultFileName = path + "/Physics/particlesDefault.xml";
  constantFileName = path + "/Physics/physConstants.xml";
  constantDefaultFileName = path + "/Physics/physDefaultConstants.xml";

  initQuantumNumberMapping();
  readFile();
}

PhysConst::~PhysConst() {
}

void PhysConst::initQuantumNumberMapping() {
  quantum_number_key_name_mapping_[QuantumNumbers::SPIN] = "spin";
  quantum_number_key_name_mapping_[QuantumNumbers::ISOSPIN] = "isospin";
  quantum_number_key_name_mapping_[QuantumNumbers::CHARGE] = "charge";
  quantum_number_key_name_mapping_[QuantumNumbers::PARITY] = "parity";
  quantum_number_key_name_mapping_[QuantumNumbers::CPARITY] = "cparity";

}

void PhysConst::readFile() {

  // Create an empty property tree object
  using boost::property_tree::ptree;
  ptree pt;
  ptree pt2;

  // Load the XML file into the property tree. If reading fails
  // (cannot open file, parse error), an exception is thrown.
  //		read_xml(particleFileName, pt);
  //first check if particles.xml existsheck if file exists
  if (FILE *file = std::fopen(particleFileName.c_str(), "r")) {
    fclose(file);
    read_xml(particleFileName, pt);
    BOOST_LOG_TRIVIAL(info)<< "PhysConst: reading particle file "<<particleFileName;
    //Otherwise try to load default file
  }
  else if (FILE *file = std::fopen(particleDefaultFileName.c_str(), "r")) {
    fclose(file);
    read_xml(particleDefaultFileName, pt);
    BOOST_LOG_TRIVIAL(info) << "PhysConst: reading particles default file "<<particleDefaultFileName;
  }
  else {
    throw std::runtime_error("Could not open default particle file!");
  }

  // Get the filename and store it in the m_file variable.
  // Note that we construct the path to the value by separating
  // the individual keys with dots. If dots appear in the keys,
  // a path type with a different separator can be used.
  // If the amplitude_setup.filename key is not found, an
  // exception is thrown.
  //    m_file = pt.get<std::string>("amplitude_setup.filename");

  // Iterate over the amplitude_setup.resonances section and
  // store all found resonance names in the m_name set.
  // The get_child() function returns a reference to the child
  // at the specified path; if there is no such child, it throws.
  // Property tree iterators are models of BidirectionalIterator.
  ParticleProperties particle_properties;
  Constant constant;

  BOOST_FOREACH( ptree::value_type const& v, pt.get_child("particleList") ) {
    if (v.first == "particle" || v.first == "particleFlatte") {
      particle_properties.id_ = v.second.get<int>("ID");
      particle_properties.name_ = v.second.get<std::string>("name");
      particle_properties.mass_ = v.second.get_child("mass").get<double>(
          "value");
      if (v.second.count("width") != 0)
        particle_properties.width_ = v.second.get_child("width").get<double>(
            "value");    //check if node "width" exists

      particle_properties.charge_ = v.second.get<int>("charge");
      if (v.second.count("isospin") != 0) {
        particle_properties.isospin_.J_numerator_ = v.second.get<unsigned int>(
            "isospin");
        particle_properties.isospin_.J_denominator_ = 1;
        particle_properties.isospin_.J_z_numerator_ = v.second.get<int>(
            "isospin_z");
      }

      particle_properties.spin_.J_numerator_ = v.second.get<unsigned int>("J");
      particle_properties.spin_.J_denominator_ = 1;

      particle_properties.parity_ = v.second.get<int>("P");
      if (v.second.count("C") != 0)
        particle_properties.cparity_ = v.second.get<int>("C");

      if (v.first == "particleFlatte") {
        //read parameters which are specific to flatte description here.

      }

      particle_properties_list_.push_back(particle_properties);

      BOOST_LOG_TRIVIAL(debug)<<"PhysConst adding particle: "<<particle_properties.name_<<" mass="<<particle_properties.mass_<<" width="<<particle_properties.width_<<" J=" <<1.0*particle_properties.spin_.J_numerator_/particle_properties.spin_.J_denominator_<<" P="<<particle_properties.parity_<< " C="<<particle_properties.cparity_;
    }
  }

//Reading XML file with physics constants
  if (FILE *file = std::fopen(constantFileName.c_str(), "r")) {
    fclose(file);
    read_xml(constantFileName, pt);
    BOOST_LOG_TRIVIAL(info)<< "PhysConst: reading file with physical constants"<<constantFileName;
    //Otherwise try to load default file
  }
  else if (FILE *file = std::fopen(constantDefaultFileName.c_str(), "r")) {
    fclose(file);
    read_xml(constantDefaultFileName, pt);
    BOOST_LOG_TRIVIAL(info) << "PhysConst: reading default file with physical constants"<<constantDefaultFileName;
  }
  else {
    throw std::runtime_error("Could not open default constants file!");
  }
//	read_xml(constantFileName, pt);
  BOOST_FOREACH( ptree::value_type const& v, pt.get_child("physConstList") ) {
    if (v.first == "constant") {
      constant.name_ = v.second.get<std::string>("name");
      constant.value_ = v.second.get_child("value").get<double>("value");
      if (v.second.count("error") != 0)
        constant.error_ = v.second.get_child("value").get<double>("error");
    }

    constants_list_.push_back(constant);

    BOOST_LOG_TRIVIAL(debug)<<"PhysConst adding constant: "<<constant.name_<<" value="<<constant.value_<<" error="<<constant.error_;
  }

  return;
}

const Constant& PhysConst::findConstant(const std::string& name) const {
  auto result = std::find_if(constants_list_.begin(), constants_list_.end(),
      [&] (const Constant& lhs) {return lhs.name_ == name;});

  if (result != constants_list_.end())
    return *result;

  std::stringstream ss;
  ss << "could not find constant with name " << name << std::endl;
  throw std::runtime_error(ss.str());
}

const ParticleProperties& PhysConst::findParticle(int pid) const {
  auto result = std::find_if(particle_properties_list_.begin(),
      particle_properties_list_.end(),
      [&] (const ParticleProperties& lhs) {return lhs.id_ == pid;});

  if (result != particle_properties_list_.end())
    return *result;

  std::stringstream ss;
  ss << "could not find particle id " << pid << std::endl;
  throw std::runtime_error(ss.str());
}

const ParticleProperties& PhysConst::findParticle(
    const std::string& name) const {
  auto result = std::find_if(particle_properties_list_.begin(),
      particle_properties_list_.end(),
      [&] (const ParticleProperties& lhs) {return lhs.name_ == name;});

  if (result != particle_properties_list_.end())
    return *result;

  std::stringstream ss;
  ss << "could not find particle name " << name << std::endl;
  throw std::runtime_error(ss.str());
}

std::vector<ParticleProperties> PhysConst::findParticlesWithQN(
    const ParticleProperties& qn) const {
  std::vector<ParticleProperties> particle_list;

  auto result = particle_properties_list_.begin();
  while (result != particle_properties_list_.end()) {
    result = std::find_if(result, particle_properties_list_.end(), qn);
    if (result != particle_properties_list_.end()) {
      particle_list.push_back(*result);
      ++result;
    }
  }

  return particle_list;
}

bool PhysConst::particleExists(const std::string& name) const {
  auto result = std::find_if(particle_properties_list_.begin(),
      particle_properties_list_.end(),
      [&] (const ParticleProperties& lhs) {return lhs.name_ == name;});

  if (result != particle_properties_list_.end())
    return true;

  return false;
}

std::string PhysConst::getQuantumNumberName(
    const QuantumNumbers& qn_type) const {
  auto result = quantum_number_key_name_mapping_.find(qn_type);
  if (result != quantum_number_key_name_mapping_.end()) {
    return result->second;
  }
  else {
    std::runtime_error("DecayGeneratorFacade::setAllowedSpinQuantumNumbers:"
        " quantum number with your specified key does not exist in the mapping!"
        " Please correct or update mapping!");
  }
}

}
