// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Core/SubSystem.hpp"

using namespace ComPWA;

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
  // Creating unique title
  std::stringstream stream;
  //  stream << "(";
  for (auto i : _recoilState)
    stream << std::to_string(i);
  stream << "_vs_";
  for (auto j = _finalStates.begin(); j != _finalStates.end(); ++j) {
    for (auto i : *j)
      stream << std::to_string(i);
    if (j != _finalStates.end() - 1)
      stream << "_";
  }

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
