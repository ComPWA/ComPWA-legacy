// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "SubSystem.hpp"

namespace ComPWA {
namespace Physics {

SubSystem::SubSystem(const boost::property_tree::ptree &pt) { load(pt); }

void SubSystem::load(const boost::property_tree::ptree &pt) {
  RecoilFinalState.clear();
  ParentRecoilFinalState.clear();
  FinalStates.clear();

  // Read subSystem definition
  auto recoil =
      pt.get_optional<std::string>("RecoilSystem.<xmlattr>.RecoilFinalState");
  if (recoil) {
    setRecoilState(stringToVectInt(recoil.get()));
  }
  auto parent_recoil = pt.get_optional<std::string>(
      "RecoilSystem.<xmlattr>.ParentRecoilFinalState");
  if (parent_recoil) {
    setParentRecoilState(stringToVectInt(parent_recoil.get()));
  }
  auto decayProducts = pt.get_child("DecayProducts");

  std::vector<std::vector<unsigned int>> tempvec;
  for (auto i : decayProducts) {
    tempvec.push_back(
        stringToVectInt(i.second.get<std::string>("<xmlattr>.FinalState")));
  }
  setFinalStates(tempvec);
}

boost::property_tree::ptree SubSystem::save() const {
  boost::property_tree::ptree pt;
  auto recoilV = getRecoilState();
  if (recoilV.size()) {
    std::string recoilStr;
    for (unsigned int i = 0; i < recoilV.size(); i++) {
      recoilStr += std::to_string(recoilV.at(i));
      if (i < recoilV.size() - 1)
        recoilStr += " ";
    }
    pt.put("RecoilSystem.<xmlattr>.FinalState", recoilStr);
  }

  boost::property_tree::ptree daughterTr;

  // Information daugher final state A
  auto finalS = getFinalStates();
  for (unsigned int j = 0; j < finalS.size(); j++) {
    std::string strA;
    if (finalS.at(j).size()) {
      for (unsigned int i = 0; i < finalS.at(j).size(); i++) {
        strA += std::to_string(finalS.at(j).at(i));
        if (i < finalS.at(j).size() - 1)
          strA += " ";
      }
    }

    boost::property_tree::ptree fTr;
    fTr.put("<xmlattr>.FinalState", strA);
    daughterTr.add_child("Particle", fTr);
  }

  pt.add_child("DecayProducts", daughterTr);

  return pt;
}

SubSystem::SubSystem(const std::vector<std::vector<unsigned int>> &FinalStates,
                     const std::vector<unsigned int> &Recoil,
                     const std::vector<unsigned int> &ParentRecoil) {
  setFinalStates(FinalStates);
  setRecoilState(Recoil);
  setParentRecoilState(ParentRecoil);
}

std::ostream &operator<<(std::ostream &stream, const SubSystem &s) {
  for (auto const &j : s.FinalStates) {
    stream << "_";
    for (auto const &i : j)
      stream << std::to_string(i);
  }
  if (s.RecoilFinalState.size() > 0) {
    stream << "_vs_";
    for (auto i : s.RecoilFinalState)
      stream << std::to_string(i);
  }

  return stream;
}

bool SubSystem::operator==(const SubSystem &b) const {
  // we assume all vectors are sorted!!!
  if (RecoilFinalState == b.RecoilFinalState && FinalStates == b.FinalStates &&
      ParentRecoilFinalState == b.ParentRecoilFinalState)
    return true;
  return false;
}

void SubSystem::setFinalStates(
    const std::vector<std::vector<unsigned int>> &v) {
  FinalStates.clear();
  for (auto fs : v) {
    std::sort(fs.begin(), fs.end());
    FinalStates.push_back(fs);
  }
}

const std::vector<std::vector<unsigned int>> &
SubSystem::getFinalStates() const {
  return FinalStates;
}

void SubSystem::setRecoilState(const std::vector<unsigned int> &r) {
  auto recoil_state_copy(r);
  std::sort(recoil_state_copy.begin(), recoil_state_copy.end());
  RecoilFinalState = recoil_state_copy;
}

void SubSystem::setParentRecoilState(const std::vector<unsigned int> &r) {
  auto recoil_state_copy(r);
  std::sort(recoil_state_copy.begin(), recoil_state_copy.end());
  ParentRecoilFinalState = recoil_state_copy;
}

const std::vector<unsigned int> &SubSystem::getRecoilState() const {
  return RecoilFinalState;
}

const std::vector<unsigned int> &SubSystem::getParentRecoilState() const {
  return ParentRecoilFinalState;
}

} // namespace Physics
} // namespace ComPWA
