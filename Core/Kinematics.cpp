// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Core/Kinematics.hpp"
#include "Core/Exceptions.hpp"
#include "Core/DataPoint.hpp"
#include "Core/Event.hpp"

using namespace ComPWA;

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
