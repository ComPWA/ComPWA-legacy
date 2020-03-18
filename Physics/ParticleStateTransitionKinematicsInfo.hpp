// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_PARTICLESTATETRANSITIONKINEMATICSINFO_HPP_
#define COMPWA_PARTICLESTATETRANSITIONKINEMATICSINFO_HPP_

#include <map>
#include <string>
#include <vector>

#include "Core/FourMomentum.hpp"
#include "Core/Properties.hpp"

namespace ComPWA {
namespace Physics {

class ParticleStateTransitionKinematicsInfo {
  std::vector<pid> InitialState;
  std::vector<pid> FinalState;
  ComPWA::ParticleList ParticleList;
  /// Four momentum of the initial particle reaction
  ComPWA::FourMomentum InitialStateP4;
  // we use a vector instead of a map here, due to cache optimizations
  std::vector<unsigned int> FinalStateEventPositionMapping;

public:
  ParticleStateTransitionKinematicsInfo(
      std::vector<pid> InitialState_, std::vector<pid> FinalState_,
      ComPWA::ParticleList ParticleList_, ComPWA::FourMomentum InitialStateP4_,
      std::vector<unsigned int> FinalStateEventPositionMapping_);

  ParticleStateTransitionKinematicsInfo(
      std::vector<pid> InitialState_, std::vector<pid> FinalState_,
      ComPWA::ParticleList ParticleList_,
      std::vector<unsigned int> FinalStateEventPositionMapping_);

  unsigned int convertFinalStateIDToPositionIndex(unsigned int fs_id) const;
  std::vector<unsigned int> convertFinalStateIDToPositionIndex(
      const std::vector<unsigned int> &fs_ids) const;

  unsigned int convertPositionIndexToFinalStateID(unsigned int pos) const;
  std::vector<unsigned int> convertPositionIndexToFinalStateID(
      const std::vector<unsigned int> &pos) const;

  double
  calculateFinalStateIDMassSum(const std::vector<unsigned int> ids) const;

  std::vector<double> getFinalStateMasses() const;
  double getInitialStateInvariantMassSquared() const;
  ComPWA::FourMomentum getInitialStateFourMomentum() const;
  unsigned int getFinalStateParticleCount() const;
  std::map<unsigned int, std::string> getFinalStateIDToNameMapping() const;
  const std::vector<pid> &getFinalStatePIDs() const { return FinalState; }

  friend std::ostream &
  operator<<(std::ostream &outstream,
             const ParticleStateTransitionKinematicsInfo &kininfo);
};

} // namespace Physics
} // namespace ComPWA

#endif
