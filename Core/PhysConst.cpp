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

using namespace boost::log;

namespace ComPWA {
PhysConst *PhysConst::_inst;

QuantumNumberTranslator::QuantumNumberTranslator() {
  quantum_number_key_name_mapping_[QuantumNumberIDs::SPIN] = "spin";
  quantum_number_key_name_mapping_[QuantumNumberIDs::ORBITAL_ANGULAR_MOMENTUM] =
      "angular-momentum";
  quantum_number_key_name_mapping_[QuantumNumberIDs::ISOSPIN] = "isospin";
  quantum_number_key_name_mapping_[QuantumNumberIDs::CHARGE] = "charge";
  quantum_number_key_name_mapping_[QuantumNumberIDs::PARITY] = "parity";
  quantum_number_key_name_mapping_[QuantumNumberIDs::CPARITY] = "cparity";

  name_quantum_number_key_mapping_["spin"] = QuantumNumberIDs::SPIN;
  name_quantum_number_key_mapping_["angular-momentum"] =
      QuantumNumberIDs::ORBITAL_ANGULAR_MOMENTUM;
  name_quantum_number_key_mapping_["isospin"] = QuantumNumberIDs::ISOSPIN;
  name_quantum_number_key_mapping_["charge"] = QuantumNumberIDs::CHARGE;
  name_quantum_number_key_mapping_["parity"] = QuantumNumberIDs::PARITY;
  name_quantum_number_key_mapping_["cparity"] = QuantumNumberIDs::CPARITY;

  quantum_number_key_type_mapping_[QuantumNumberIDs::SPIN] =
      QuantumNumberType::SPIN_LIKE;
  quantum_number_key_type_mapping_[QuantumNumberIDs::ORBITAL_ANGULAR_MOMENTUM] =
      QuantumNumberType::SPIN_LIKE;
  quantum_number_key_type_mapping_[QuantumNumberIDs::ISOSPIN] =
      QuantumNumberType::SPIN_LIKE;
  quantum_number_key_type_mapping_[QuantumNumberIDs::CHARGE] =
      QuantumNumberType::INTEGER_LIKE;
  quantum_number_key_type_mapping_[QuantumNumberIDs::PARITY] =
      QuantumNumberType::INTEGER_LIKE;
  quantum_number_key_type_mapping_[QuantumNumberIDs::CPARITY] =
      QuantumNumberType::INTEGER_LIKE;
};

QuantumNumberTranslator::~QuantumNumberTranslator(){};

QuantumNumberType QuantumNumberTranslator::getQuantumNumberType(
    const std::string &qn_name) const {
  auto result =
      quantum_number_key_type_mapping_.find(getQuantumNumberEnum(qn_name));
  if (result != quantum_number_key_type_mapping_.end()) {
    return result->second;
  } else {
    std::runtime_error(
        "QuantumNumberTranslator::getQuantumNumberType:"
        " quantum number with your specified key does not exist in the mapping!"
        " Please correct or update mapping!");
  }
}

std::string QuantumNumberTranslator::getQuantumNumberName(
    const QuantumNumberIDs &qn_type) const {
  auto result = quantum_number_key_name_mapping_.find(qn_type);
  if (result != quantum_number_key_name_mapping_.end()) {
    return result->second;
  } else {
    std::runtime_error(
        "QuantumNumberTranslator::getQuantumNumberName:"
        " quantum number with your specified key does not exist in the mapping!"
        " Please correct or update mapping!");
  }
}

QuantumNumberIDs QuantumNumberTranslator::getQuantumNumberEnum(
    const std::string &qn_name) const {
  auto result = name_quantum_number_key_mapping_.find(qn_name);
  if (result != name_quantum_number_key_mapping_.end()) {
    return result->second;
  } else {
    std::runtime_error(
        "QuantumNumberTranslator::getQuantumNumberEnum:"
        " quantum number with your specified key does not exist in the mapping!"
        " Please correct or update mapping!");
  }
}

ComPWA::Spin SpinWave::getSpinLikeQuantumNumber(QuantumNumberIDs qn_id) const {
  auto spin_result = spin_like_quantum_numbers_.find(
      QuantumNumberTranslator::Instance().getQuantumNumberName(qn_id));
  if (spin_result != spin_like_quantum_numbers_.end())
    return spin_result->second;
  else {
    std::stringstream ss;
    ss << "SpinWave::getSpinLikeQuantumNumber: did not find quantum number ";
    ss << QuantumNumberTranslator::Instance().getQuantumNumberName(qn_id);
    ss << " in spin like qn list!";
    std::runtime_error(ss.str());
  }
}

int SpinWave::getIntLikeQuantumNumber(QuantumNumberIDs qn_id) const {
  auto int_result = integer_like_quantum_numbers_.find(
      QuantumNumberTranslator::Instance().getQuantumNumberName(qn_id));
  if (int_result != integer_like_quantum_numbers_.end())
    return int_result->second;
  else {
    std::stringstream ss;
    ss << "SpinWave::getIntLikeQuantumNumber: did not find quantum number ";
    ss << QuantumNumberTranslator::Instance().getQuantumNumberName(qn_id);
    ss << " in int like qn list!";
    std::runtime_error(ss.str());
  }
}

double SpinWave::getDoubleLikeQuantumNumber(QuantumNumberIDs qn_id) const {
  auto double_result = double_like_quantum_numbers_.find(
      QuantumNumberTranslator::Instance().getQuantumNumberName(qn_id));
  if (double_result != double_like_quantum_numbers_.end())
    return double_result->second;
  else {
    std::stringstream ss;
    ss << "SpinWave::getDoubleLikeQuantumNumber: did not find quantum number ";
    ss << QuantumNumberTranslator::Instance().getQuantumNumberName(qn_id);
    ss << " in double like qn list!";
    std::runtime_error(ss.str());
  }
}

PhysConst::PhysConst(std::string path) {
  particleFileName = path + "particles.xml";
  particleDefaultFileName = path + "particlesDefault.xml";
  constantFileName = path + "physConstants.xml";
  constantDefaultFileName = path + "physDefaultConstants.xml";

  readFile();
}

PhysConst::~PhysConst() {}

void PhysConst::readFile() {
  // Create an empty property tree object
  using boost::property_tree::ptree;
  ptree pt;
  ptree pt2;

  // Load the XML file into the property tree. If reading fails
  // (cannot open file, parse error), an exception is thrown.
  //		read_xml(particleFileName, pt);
  // first check if particles.xml existsheck if file exists
  if (FILE *file = std::fopen(particleFileName.c_str(), "r")) {
    fclose(file);
    read_xml(particleFileName, pt);
    LOG(info) << "PhysConst: reading particle file " << particleFileName;
    // Otherwise try to load default file
  } else if (FILE *file = std::fopen(particleDefaultFileName.c_str(), "r")) {
    fclose(file);
    read_xml(particleDefaultFileName, pt);
    LOG(info) << "PhysConst: reading particles default file "
              << particleDefaultFileName;
  } else {
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

  for (auto const &v : pt.get_child("particleList")) {
    if (v.first == "particle" || v.first == "particleFlatte") {
      particle_properties.id_ = v.second.get<int>("ID");
      particle_properties.name_ = v.second.get<std::string>("name");
      particle_properties.mass_ =
          v.second.get_child("mass").get<double>("value");
      if (v.second.count("width") != 0)
        particle_properties.width_ = v.second.get_child("width").get<double>(
            "value"); // check if node "width" exists

      for (auto const &qn : v.second.get_child("quantum_numbers")) {
        QuantumNumberType qn_type =
            QuantumNumberTranslator::Instance().getQuantumNumberType(qn.first);
        if (qn_type == QuantumNumberType::SPIN_LIKE) {
          ComPWA::Spin s;
          if (QuantumNumberTranslator::Instance().getQuantumNumberEnum(
                  qn.first) == QuantumNumberIDs::SPIN)
            s.SetUseZ(false);
          s.SetNumerator(qn.second.get<unsigned int>("numerator"));
          s.SetDenominator(qn.second.get<unsigned int>("denominator"));
          if (qn.second.count("z_numerator") != 0)
            s.SetZNumerator(qn.second.get<int>("z_numerator"));
          particle_properties.spin_like_quantum_numbers_[qn.first] = s;
        } else if (qn_type == QuantumNumberType::INTEGER_LIKE) {
          particle_properties.integer_like_quantum_numbers_[qn.first] =
              v.second.get_child("quantum_numbers").get<int>(qn.first);
        } else if (qn_type == QuantumNumberType::DOUBLE_LIKE) {
          particle_properties.double_like_quantum_numbers_[qn.first] =
              v.second.get_child("quantum_numbers").get<double>(qn.first);
        }
      }

      if (v.first == "particleFlatte") {
        // read parameters which are specific to flatte description here.
      }

      particle_properties_list_.push_back(particle_properties);

      ComPWA::Spin s =
          particle_properties.getSpinLikeQuantumNumber(QuantumNumberIDs::SPIN);
      LOG(debug) << "PhysConst adding particle: " << particle_properties.name_
                 << " mass=" << particle_properties.mass_
                 << " width=" << particle_properties.width_
                 << " J=" << s.GetSpin() << " P="
                 << particle_properties.getIntLikeQuantumNumber(
                        QuantumNumberIDs::PARITY)
                 << " C="
                 << particle_properties.getIntLikeQuantumNumber(
                        QuantumNumberIDs::CPARITY);
    }
  }

  // Reading XML file with physics constants
  if (FILE *file = std::fopen(constantFileName.c_str(), "r")) {
    fclose(file);
    read_xml(constantFileName, pt);
    LOG(info) << "PhysConst: reading file with physical constants"
              << constantFileName;
    // Otherwise try to load default file
  } else if (FILE *file = std::fopen(constantDefaultFileName.c_str(), "r")) {
    fclose(file);
    read_xml(constantDefaultFileName, pt);
    LOG(info) << "PhysConst: reading default file with physical constants"
              << constantDefaultFileName;
  } else {
    throw std::runtime_error("Could not open default constants file!");
  }
  //	read_xml(constantFileName, pt);
  BOOST_FOREACH (ptree::value_type const &v, pt.get_child("physConstList")) {
    if (v.first == "constant") {
      constant.name_ = v.second.get<std::string>("name");
      constant.value_ = v.second.get_child("value").get<double>("value");
      if (v.second.count("error") != 0)
        constant.error_ = v.second.get_child("value").get<double>("error");
    }

    constants_list_.push_back(constant);

    LOG(debug) << "PhysConst adding constant: " << constant.name_
               << " value=" << constant.value_ << " error=" << constant.error_;
  }

  return;
}

const Constant &PhysConst::findConstant(const std::string &name) const {
  auto result =
      std::find_if(constants_list_.begin(), constants_list_.end(),
                   [&](const Constant &lhs) { return lhs.name_ == name; });

  if (result != constants_list_.end())
    return *result;

  std::stringstream ss;
  ss << "could not find constant with name " << name << std::endl;
  throw std::runtime_error(ss.str());
}

const ParticleProperties &PhysConst::findParticle(int pid) const {
  auto result = std::find_if(
      particle_properties_list_.begin(), particle_properties_list_.end(),
      [&](const ParticleProperties &lhs) { return lhs.id_ == pid; });

  if (result != particle_properties_list_.end())
    return *result;

  std::stringstream ss;
  ss << "could not find particle id " << pid << std::endl;
  throw std::runtime_error(ss.str());
}

const ParticleProperties &
PhysConst::findParticle(const std::string &name) const {
  auto result = std::find_if(
      particle_properties_list_.begin(), particle_properties_list_.end(),
      [&](const ParticleProperties &lhs) { return lhs.name_ == name; });

  if (result != particle_properties_list_.end())
    return *result;

  std::stringstream ss;
  ss << "could not find particle name " << name << std::endl;
  throw std::runtime_error(ss.str());
}

std::vector<ParticleProperties>
PhysConst::findParticlesWithQN(const ParticleProperties &qn) const {
  std::vector<ParticleProperties> particle_list;

  auto result = particle_properties_list_.begin();
  while (result != particle_properties_list_.end()) {
    result = std::find(result, particle_properties_list_.end(), qn);
    if (result != particle_properties_list_.end()) {
      particle_list.push_back(*result);
      ++result;
    }
  }

  return particle_list;
}

bool PhysConst::particleExists(const std::string &name) const {
  auto result = std::find_if(
      particle_properties_list_.begin(), particle_properties_list_.end(),
      [&](const ParticleProperties &lhs) { return lhs.name_ == name; });

  if (result != particle_properties_list_.end())
    return true;

  return false;
}

} /* namespace ComPWA */
