// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Core/SubSystem.hpp"

using namespace ComPWA;

SubSystem::SubSystem(const boost::property_tree::ptree &pt) { load(pt); }

void SubSystem::load(const boost::property_tree::ptree &pt) {
  _recoilState.clear();
  _finalStates.clear();
  helicities_.clear();
  _finalStatesNames.clear();
  
  // Read subSystem definition
  auto recoil =
      pt.get_optional<std::string>("RecoilSystem.<xmlattr>.FinalState");
  if (recoil) {
    _recoilState = stringToVectInt(recoil.get());
  }
  auto decayProducts = pt.get_child("DecayProducts");

  for (auto i : decayProducts) {
    _finalStates.push_back(
        stringToVectInt(i.second.get<std::string>("<xmlattr>.FinalState")));
    _finalStatesNames.push_back(i.second.get<std::string>("<xmlattr>.Name"));
    helicities_.push_back(i.second.get<int>("<xmlattr>.Helicity"));
  }
}

boost::property_tree::ptree SubSystem::save() const {
  boost::property_tree::ptree pt;
  auto recoilV = GetRecoilState();
  if (recoilV.size()) {
    std::string recoilStr;
    for (auto i = 0; i < recoilV.size(); i++) {
      recoilStr += std::to_string(recoilV.at(i));
      if (i < recoilV.size() - 1)
        recoilStr += " ";
    }
    pt.put("RecoilSystem.<xmlattr>.FinalState", recoilStr);
  }

  boost::property_tree::ptree daughterTr;

  // Information daugher final state A
  auto finalS = GetFinalStates();
  auto helicities = GetHelicities();
  auto finalSNames = GetFinalStatesNames();
  for (int j = 0; j < finalS.size(); j++) {
    std::string strA;
    if (finalS.at(j).size()) {
      for (auto i = 0; i < finalS.at(j).size(); i++) {
        strA += std::to_string(finalS.at(j).at(i));
        if (i < finalS.at(j).size() - 1)
          strA += " ";
      }
    }

    boost::property_tree::ptree fTr;
    fTr.put("<xmlattr>.Name", finalSNames.at(j));
    fTr.put("<xmlattr>.FinalState", strA);
    fTr.put("<xmlattr>.Helicity", helicities.at(j));
    daughterTr.add_child("Particle", fTr);
  }

  pt.add_child("DecayProducts", daughterTr);

  return pt;
}

SubSystem::SubSystem(std::vector<int> recoilS,
                     std::vector<std::vector<int>> finalStates)
    : _recoilState(recoilS), _finalStates(finalStates) {

  _title = to_string();
  // LOG(trace) << "SubSystem::SubSystem() | Creating sub system "<<title;
}

/// Constructor for an isobar model in which it is assumed that particles
/// always decay to two final states.
SubSystem::SubSystem(std::vector<int> recoilS, std::vector<int> finalA,
                     std::vector<int> finalB)
    : _recoilState(recoilS) {
  std::vector<std::vector<int>> tmp;
  tmp.push_back(finalA);
  tmp.push_back(finalB);
  SetFinalStates(tmp);

  _title = to_string();
  // LOG(trace) << "SubSystem::SubSystem() | Creating sub system "<<title;
}

std::string SubSystem::to_string() const {
  std::stringstream stream;

  for (auto j = _finalStates.begin(); j != _finalStates.end(); ++j) {
    for (auto i : *j)
      stream << std::to_string(i);
    if (j != _finalStates.end() - 1)
      stream << "_";
  }
  stream << "_vs_";
  for (auto i : _recoilState)
    stream << std::to_string(i);

  return stream.str();
}

bool SubSystem::operator==(const SubSystem &b) const {
  if (_recoilState == b._recoilState && _finalStates == b._finalStates)
    return true;
  return false;
}

void SubSystem::SetFinalStates(std::vector<std::vector<int>> v) {
  _finalStates = v;
}

const std::vector<std::string> &SubSystem::GetFinalStatesNames() const {
  return _finalStatesNames;
}

void SubSystem::SetFinalStatesNames(std::vector<std::string> n) {
  if (n.size() != _finalStates.size()) {
    throw std::runtime_error("SubSystem::SetFinalStatesNames() | Length of "
                             "vectors does not match with the number of "
                             "final states.");
  }
  _finalStatesNames = n;
}

const std::vector<std::vector<int>> &SubSystem::GetFinalStates() const {
  return _finalStates;
}

const std::vector<int> SubSystem::GetHelicities() const {
  if (helicities_.size() != _finalStates.size())
    throw std::runtime_error("SubSystem::GetHelicities() | Helicities are "
                             "not defined for all final states!");
  return helicities_;
}
void SubSystem::SetHelicities(std::vector<int> hel) {
  if (hel.size() != _finalStates.size())
    throw std::runtime_error("SubSystem::SetHelicities() | Helicities are "
                             "not defined for all final states!");
  helicities_ = hel;
}

void SubSystem::SetRecoilState(const std::vector<int> r) { _recoilState = r; }

const std::vector<int> &SubSystem::GetRecoilState() const {
  return _recoilState;
}
