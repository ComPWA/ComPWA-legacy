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

PhysConst::PhysConst(std::string filename) {

  _constList.push_back(Constant("Pi", 3.141592654));
  _constList.push_back(Constant("c", 299792458)); /* Units: m/s */

  // Create an empty property tree object
  boost::property_tree::ptree pt;
  try {
    read_xml(_particleFileName, pt);
    LOG(info) << "PhysConst: reading particle file " << _particleFileName;
  } catch (std::exception &ex) {
    throw;
  }
  readTree(pt);
}

PhysConst::PhysConst(boost::property_tree::ptree pt) {

  _constList.push_back(Constant("Pi", 3.141592654));
  _constList.push_back(Constant("c", 299792458)); /* Units: m/s */

  readTree(pt);
}

PhysConst::~PhysConst() {}

void PhysConst::readTree(boost::property_tree::ptree pt) {

  for (auto const &v : pt.get_child("ParticleList")) {
    _partList.push_back(ParticleProperties(v.second));
    auto last = _partList.back();
    LOG(debug) << "PhysConst adding particle: " << last.GetName()
               << " mass=" << last.GetMass() << " J=" << last.GetSpin()
               << " P=" << last.GetParity() << " C=" << last.GetCparity();
  }

  for (auto const &v : pt.get_child("PhysConstantsList")) {
    _constList.push_back(Constant(v.second));
    auto last = _constList.back();
    LOG(debug) << "PhysConst adding constant: " << last.GetName()
               << " value=" << last.GetValue();
  }

  return;
}

const ComPWA::Constant &PhysConst::FindConstant(const std::string &name) const {
  auto result =
      std::find_if(_constList.begin(), _constList.end(),
                   [&](const Constant &lhs) { return lhs.GetName() == name; });

  if (result == _constList.end()) {
    std::stringstream ss;
    ss << "PhysConst::findConstant() | Constant " << name
       << " not found in list!";
    throw std::runtime_error(ss.str());
  }
  return *result;
}

const ParticleProperties &PhysConst::FindParticle(int pid) const {
  auto result = std::find_if(
      _partList.begin(), _partList.end(),
      [&](const ParticleProperties &lhs) { return lhs.GetId() == pid; });

  if (result == _partList.end()) {
    std::stringstream ss;
    ss << "PhysConst::findParticle() | Particle " << pid
       << " not found in list!";
    throw std::runtime_error(ss.str());
  }
  return *result;
}

const ParticleProperties &
PhysConst::FindParticle(const std::string &name) const {
  auto result = std::find_if(
      _partList.begin(), _partList.end(),
      [&](const ParticleProperties &lhs) { return lhs.GetName() == name; });

  if (result == _partList.end()) {
    std::stringstream ss;
    ss << "PhysConst::findParticle() | Particle " << name
       << " not found in list!";
    throw std::runtime_error(ss.str());
  }
  return *result;
}

// std::vector<ParticleProperties>
// PhysConst::findParticlesWithQN(const ParticleProperties &qn) const {
//  std::vector<ParticleProperties> particle_list;
//
//  auto result = _partList.begin();
//  while (result != _partList.end()) {
//    result = std::find(result, _partList.end(), qn);
//    if (result != _partList.end()) {
//      particle_list.push_back(*result);
//      ++result;
//    }
//  }
//
//  return particle_list;
//}

bool PhysConst::ParticleExists(const std::string &name) const {
  try {
    FindParticle(name);
  } catch (std::exception &ex) {
    return false;
  }
  return true;
}

} /* namespace ComPWA */
