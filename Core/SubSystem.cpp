// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Core/SubSystem.hpp"

using namespace ComPWA;

SubSystem::SubSystem(const boost::property_tree::ptree &pt) { load(pt); }

void SubSystem::load(const boost::property_tree::ptree &pt) {
  RecoilFinalState.clear();
  ParentRecoilFinalState.clear();
  FinalStates.clear();
  Helicities.clear();
  FinalStateNames.clear();

  // Read subSystem definition
  auto recoil =
      pt.get_optional<std::string>("RecoilSystem.<xmlattr>.RecoilFinalState");
  if (recoil) {
    RecoilFinalState = stringToVectInt(recoil.get());
    std::sort(RecoilFinalState.begin(), RecoilFinalState.end());
  }
  auto parent_recoil = pt.get_optional<std::string>(
      "RecoilSystem.<xmlattr>.ParentRecoilFinalState");
  if (parent_recoil) {
    ParentRecoilFinalState = stringToVectInt(parent_recoil.get());
    std::sort(ParentRecoilFinalState.begin(), ParentRecoilFinalState.end());
  }
  auto decayProducts = pt.get_child("DecayProducts");

  for (auto i : decayProducts) {
    auto tempvec(
        stringToVectInt(i.second.get<std::string>("<xmlattr>.FinalState")));
    std::sort(tempvec.begin(), tempvec.end());
    FinalStates.push_back(
        stringToVectInt(i.second.get<std::string>("<xmlattr>.FinalState")));
    FinalStateNames.push_back(i.second.get<std::string>("<xmlattr>.Name"));
    Helicities.push_back(i.second.get<double>("<xmlattr>.Helicity"));
  }
  // std::sort(_finalStates.begin(), _finalStates.end());
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
  auto helicities = getHelicities();
  auto finalSNames = getFinalStatesNames();
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
    fTr.put("<xmlattr>.Name", finalSNames.at(j));
    fTr.put("<xmlattr>.FinalState", strA);
    fTr.put("<xmlattr>.Helicity", helicities.at(j));
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
  Title = to_string();
  // LOG(TRACE) << "SubSystem::SubSystem() | Creating sub system "<<title;
}

std::string SubSystem::to_string() const {
  std::stringstream stream;

  for (auto const &j : FinalStates) {
    for (auto const &i : j)
      stream << std::to_string(i);
    stream << "_";
  }
  stream << "vs_";
  for (auto i : RecoilFinalState)
    stream << std::to_string(i);

  return stream.str();
}

bool SubSystem::operator==(const SubSystem &b) const {
  // we assume all vectors are sorted!!!
  if (RecoilFinalState == b.RecoilFinalState && FinalStates == b.FinalStates)
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
  // std::sort(_finalStates.begin(), _finalStates.end());
}

const std::vector<std::string> &SubSystem::getFinalStatesNames() const {
  return FinalStateNames;
}

void SubSystem::setFinalStatesNames(const std::vector<std::string> &n) {
  if (n.size() != FinalStates.size()) {
    throw std::runtime_error("SubSystem::SetFinalStatesNames() | Length of "
                             "vectors does not match with the number of "
                             "final states.");
  }
  FinalStateNames = n;
}

const std::vector<std::vector<unsigned int>> &
SubSystem::getFinalStates() const {
  return FinalStates;
}

const std::vector<double> SubSystem::getHelicities() const {
  if (Helicities.size() != FinalStates.size())
    throw std::runtime_error("SubSystem::GetHelicities() | Helicities are "
                             "not defined for all final states!");
  return Helicities;
}
void SubSystem::setHelicities(const std::vector<double> &hel) {
  if (hel.size() != FinalStates.size())
    throw std::runtime_error("SubSystem::SetHelicities() | Helicities are "
                             "not defined for all final states!");
  Helicities = hel;
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
