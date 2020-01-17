// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "ParticleStateTransitionKinematicsInfo.hpp"

namespace ComPWA {
namespace Physics {

ParticleStateTransitionKinematicsInfo::ParticleStateTransitionKinematicsInfo(
    std::vector<pid> InitialState_, std::vector<pid> FinalState_,
    ComPWA::ParticleList ParticleList_, ComPWA::FourMomentum InitialStateP4_,
    std::vector<unsigned int> FinalStateEventPositionMapping_)
    : InitialState(InitialState_), FinalState(FinalState_),
      ParticleList(ParticleList_), InitialStateP4(InitialStateP4_),
      FinalStateEventPositionMapping(FinalStateEventPositionMapping_) {
  // If the cms four-momentum is not set of set it here
  if (InitialStateP4 == FourMomentum(0, 0, 0, 0) && InitialState.size() == 1) {
    double sqrtS =
        findParticle(ParticleList, InitialState.at(0)).getMass().Value;
    InitialStateP4 = ComPWA::FourMomentum(0, 0, 0, sqrtS);
  } // Make sure cms momentum is set
  else if (InitialStateP4 == FourMomentum(0, 0, 0, 0))
    assert(false);
}

ParticleStateTransitionKinematicsInfo::ParticleStateTransitionKinematicsInfo(
    std::vector<pid> InitialState_, std::vector<pid> FinalState_,
    ComPWA::ParticleList ParticleList_,
    std::vector<unsigned int> FinalStateEventPositionMapping_)
    : ParticleStateTransitionKinematicsInfo(
          InitialState_, FinalState_, ParticleList_,
          [&]() {
            if (InitialState_.size() == 1) {
              double sqrtS = findParticle(ParticleList_, InitialState_.at(0))
                                 .getMass()
                                 .Value;
              return ComPWA::FourMomentum(0, 0, 0, sqrtS);
            } // Make sure cms momentum is set
            else {
              throw std::runtime_error(
                  "ParticleStateTransitionKinematicsInfo(): constructing the "
                  "info without a initial state four momentum is only possible "
                  "with a single particle initial state!");
            }
          }(),
          FinalStateEventPositionMapping_) {}

unsigned int
ParticleStateTransitionKinematicsInfo::convertFinalStateIDToPositionIndex(
    unsigned int fs_id) const {
  const auto &fsepMapping(FinalStateEventPositionMapping);
  auto result = std::find(fsepMapping.begin(), fsepMapping.end(), fs_id);
  return std::distance(fsepMapping.begin(), result);
}

std::vector<unsigned int>
ParticleStateTransitionKinematicsInfo::convertFinalStateIDToPositionIndex(
    const std::vector<unsigned int> &fs_ids) const {
  std::vector<unsigned int> pos_indices;
  pos_indices.reserve(fs_ids.size());
  for (auto fs_id : fs_ids) {
    pos_indices.push_back(convertFinalStateIDToPositionIndex(fs_id));
  }
  return pos_indices;
}

unsigned int
ParticleStateTransitionKinematicsInfo::convertPositionIndexToFinalStateID(
    unsigned int pos) const {
  return FinalStateEventPositionMapping[pos];
}

std::vector<unsigned int>
ParticleStateTransitionKinematicsInfo::convertPositionIndexToFinalStateID(
    const std::vector<unsigned int> &pos) const {
  std::vector<unsigned int> fsids;
  for (auto x : pos)
    fsids.push_back(convertPositionIndexToFinalStateID(x));
  return fsids;
}

double ParticleStateTransitionKinematicsInfo::calculateFinalStateIDMassSum(
    const std::vector<unsigned int> ids) const {
  double MassSum(0.0);
  for (auto i : ids) {
    unsigned int index = convertFinalStateIDToPositionIndex(i);
    MassSum += findParticle(ParticleList, FinalState.at(index)).getMass().Value;
  }
  return MassSum;
}

std::vector<double>
ParticleStateTransitionKinematicsInfo::getFinalStateMasses() const {
  std::vector<double> FinalStateMasses;
  for (auto ParticlePid : FinalState) { // particle 0 is mother particle
    FinalStateMasses.push_back(
        findParticle(ParticleList, ParticlePid).getMass().Value);
  }
  return FinalStateMasses;
}

double
ParticleStateTransitionKinematicsInfo::getInitialStateInvariantMassSquared()
    const {
  return InitialStateP4.invariantMassSquared();
}

ComPWA::FourMomentum
ParticleStateTransitionKinematicsInfo::getInitialStateFourMomentum() const {
  return InitialStateP4;
}

unsigned int
ParticleStateTransitionKinematicsInfo::getFinalStateParticleCount() const {
  return FinalState.size();
}

std::map<unsigned int, std::string>
ParticleStateTransitionKinematicsInfo::getFinalStateIDToNameMapping() const {
  std::map<unsigned int, std::string> mapping;
  for (unsigned int i = 0; i < FinalState.size(); ++i) {
    mapping[FinalStateEventPositionMapping[i]] =
        findParticle(ParticleList, FinalState[i]).getName();
  }
  return mapping;
}

std::ostream &operator<<(std::ostream &outstream,
                         const ParticleStateTransitionKinematicsInfo &kininfo) {
  // Create title
  outstream << "( ";
  for (auto i : kininfo.InitialState)
    outstream << findParticle(kininfo.ParticleList, i).getName() << " ";
  outstream << ")->( ";
  for (unsigned int i = 0; i < kininfo.FinalState.size(); ++i)
    outstream
        << findParticle(kininfo.ParticleList, kininfo.FinalState[i]).getName()
        << "[ID=" << kininfo.FinalStateEventPositionMapping[i] << "] ";
  outstream << ")";

  outstream << "\nEvent position to final state ID mapping:\n";
  for (unsigned int i = 0; i < kininfo.FinalStateEventPositionMapping.size();
       ++i) {
    outstream << i << ": " << kininfo.FinalStateEventPositionMapping[i] << "\n";
  }

  return outstream;
}

} // namespace Physics
} // namespace ComPWA
