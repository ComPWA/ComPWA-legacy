// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_PARTICLESTATETRANSITIONKINEMATICSINFO_HPP_
#define COMPWA_PARTICLESTATETRANSITIONKINEMATICSINFO_HPP_

#include <map>
#include <string>
#include <vector>

#include "Core/Particle.hpp"
#include "Core/Properties.hpp"

namespace ComPWA {
namespace Physics {

class ParticleStateTransitionKinematicsInfo {
  std::vector<pid> InitialState;
  std::vector<pid> FinalState;
  std::shared_ptr<PartList> ParticleList;
  /// Four momentum of the initial particle reaction
  ComPWA::FourMomentum InitialStateP4;
  // we use a vector instead of a map here, due to cache optimizations
  std::vector<unsigned int> FinalStateEventPositionMapping;

public:
  ParticleStateTransitionKinematicsInfo(
      const std::vector<pid> &InitialState_,
      const std::vector<pid> &FinalState_,
      std::shared_ptr<PartList> ParticleList_,
      const ComPWA::FourMomentum &InitialStateP4_,
      const std::vector<unsigned int> &FinalStateEventPositionMapping_);

  ParticleStateTransitionKinematicsInfo(
      const std::vector<pid> &InitialState_,
      const std::vector<pid> &FinalState_,
      std::shared_ptr<PartList> ParticleList_,
      const std::vector<unsigned int> &FinalStateEventPositionMapping_);

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

  friend std::ostream &
  operator<<(std::ostream &outstream,
             const ParticleStateTransitionKinematicsInfo &kininfo);
};

} // namespace Physics
} // namespace ComPWA

#endif
