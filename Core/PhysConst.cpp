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

PhysConst::PhysConst(std::string filename) {

  _constList.push_back(Constant("Pi", M_PI));
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

  auto particleTree = pt.get_child_optional("ParticleList");
  if (particleTree) {
    for (auto const &v : particleTree.get()) {
      _partList.push_back(ParticleProperties(v.second));
      auto last = _partList.back();
      LOG(debug) << "PhysConst::readTree() | Adding particle " << last.GetName()
                 << " (id=" << last.GetId() << ") "
                 << " J(PC)=" << last.GetSpin() << "(" << last.GetParity()
                 << last.GetCparity() << ") "
                 << " mass=" << last.GetMass()
                 << " decayType=" << last.GetDecayType();
    }
  }

  auto constTree = pt.get_child_optional("PhysConstantsList");
  if (constTree) {
    for (auto const &v : constTree.get()) {
      _constList.push_back(Constant(v.second));
      auto last = _constList.back();
      LOG(debug) << "PhysConst::readTree() | Adding constant " << last.GetName()
                 << " value=" << last.GetValue();
    }
  }

  return;
}

const ComPWA::Constant &PhysConst::FindConstant(const std::string &name) const {
  auto result =
      std::find_if(_constList.begin(), _constList.end(),
                   [&](const Constant &lhs) { return lhs.GetName() == name; });

  if (result == _constList.end()) {
    std::stringstream ss;
    ss << "PhysConst::FindConstant() | Constant " << name
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
    ss << "PhysConst::FindParticle() | Particle id=" << pid
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
    ss << "PhysConst::FindParticle() | Particle " << name
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
