// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "SubSystem.hpp"

namespace ComPWA {
namespace Physics {

SubSystem::SubSystem(const std::vector<std::vector<unsigned int>> &FinalStates,
                     const std::vector<unsigned int> &Recoil,
                     const std::vector<unsigned int> &ParentRecoil)
    : DecayProductFinalStates(FinalStates), RecoilFinalState(Recoil),
      ParentRecoilFinalState(ParentRecoil) {

  for (auto &fs : DecayProductFinalStates) {
    std::sort(fs.begin(), fs.end());
  }
  std::sort(RecoilFinalState.begin(), RecoilFinalState.end());
  std::sort(ParentRecoilFinalState.begin(), ParentRecoilFinalState.end());
}

std::ostream &operator<<(std::ostream &stream, const SubSystem &s) {
  for (auto const &j : s.DecayProductFinalStates) {
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
  if (RecoilFinalState == b.RecoilFinalState &&
      DecayProductFinalStates == b.DecayProductFinalStates &&
      ParentRecoilFinalState == b.ParentRecoilFinalState)
    return true;
  return false;
}

const std::vector<std::vector<unsigned int>> &
SubSystem::getFinalStates() const {
  return DecayProductFinalStates;
}

const std::vector<unsigned int> &SubSystem::getRecoilState() const {
  return RecoilFinalState;
}

const std::vector<unsigned int> &SubSystem::getParentRecoilState() const {
  return ParentRecoilFinalState;
}

} // namespace Physics
} // namespace ComPWA
