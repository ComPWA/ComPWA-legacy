// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

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

PhysConst *PhysConst::_inst = 0;

PhysConst::PhysConst(std::string filename) {

  _constList.push_back(Constant("Pi", M_PI));
  _constList.push_back(Constant("c", 299792458)); // Units: m/s

  // Create an empty property tree object
  boost::property_tree::ptree pt;
  try {
    read_xml(filename, pt);
    LOG(info) << "PhysConst: reading particle file " << filename;
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
      
      //cparity is optional
      double cparity = 0.0;
      try { cparity = last.GetQuantumNumber("Cparity"); }
      catch (std::exception& ex) { }
      
      LOG(info) << "PhysConst::readTree() | Adding particle " << last.GetName()
                 << " (id=" << last.GetId() << ") "
                 << " J(PC)=" << last.GetSpinQuantumNumber("Spin")
                 << "(" << last.GetQuantumNumber("Parity")
                 <<  cparity<< ") "
                 << " mass=" << last.GetMass()
                 << " decayType=" << last.GetDecayType();
    }
  }

  auto constTree = pt.get_child_optional("PhysConstantsList");
  if (constTree) {
    for (auto const &v : constTree.get()) {
      _constList.push_back(Constant(v.second));
      auto last = _constList.back();
      LOG(info) << "PhysConst::readTree() | Adding constant " << last.GetName()
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

const ParticleProperties &PhysConst::FindParticle(pid id) const {
  auto result = std::find_if(
      _partList.begin(), _partList.end(),
      [&](const ParticleProperties &lhs) { return lhs.GetId() == id; });

  if (result == _partList.end()) {
    std::stringstream ss;
    ss << "PhysConst::FindParticle() | Particle id=" << id
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

bool PhysConst::ParticleExists(const std::string &name) const {
  try {
    FindParticle(name);
  } catch (std::exception &ex) {
    return false;
  }
  return true;
}

} /* namespace ComPWA */
