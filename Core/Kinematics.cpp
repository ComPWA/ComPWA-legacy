// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Core/Kinematics.hpp"

namespace ComPWA {

Kinematics::Kinematics(
    const ComPWA::KinematicsProperties &KinematicsProperties_)
    : KinematicsProperties(KinematicsProperties_), HasPhspVolume(false),
      PhspVolume(1.0) {
  // If the cms four-momentum is not set of set it here
  if (KinematicsProperties.InitialStateP4 == FourMomentum(0, 0, 0, 0) &&
      KinematicsProperties.InitialState.size() == 1) {
    double sqrtS = FindParticle(KinematicsProperties.ParticleList,
                                KinematicsProperties.InitialState.at(0))
                       .GetMass();
    KinematicsProperties.InitialStateP4 = ComPWA::FourMomentum(0, 0, 0, sqrtS);
  }
  // Make sure cms momentum is set
  if (KinematicsProperties.InitialStateP4 == FourMomentum(0, 0, 0, 0))
    assert(false);
}

double Kinematics::phspVolume() const {
  if (!HasPhspVolume) {
    const_cast<double &>(PhspVolume) = calculatePhspVolume();
    const_cast<bool &>(HasPhspVolume) = true;
  }
  return PhspVolume;
}

void Kinematics::setPhspVolume(double vol) {
  PhspVolume = vol;
  HasPhspVolume = true;
  
  LOG(INFO)<<"Kinematics::setPhspVolume() | Setting phase space "
  "volume to "<<std::to_string(phspVolume())<<".";
}

unsigned int Kinematics::convertFinalStateIDToPositionIndex(
    unsigned int fs_id) const {
  const auto &fsepMapping(KinematicsProperties.FinalStateEventPositionMapping);
  auto result = std::find(fsepMapping.begin(), fsepMapping.end(), fs_id);
  return std::distance(fsepMapping.begin(), result);
}

std::vector<unsigned int>
Kinematics::convertFinalStateIDToPositionIndex(
    const std::vector<unsigned int> &fs_ids) const {
  std::vector<unsigned int> pos_indices;
  pos_indices.reserve(fs_ids.size());
  for (auto fs_id : fs_ids) {
    pos_indices.push_back(convertFinalStateIDToPositionIndex(fs_id));
  }
  return pos_indices;
}


}
