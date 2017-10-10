// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Core/Kinematics.hpp"
#include "Core/Exceptions.hpp"
#include "Core/DataPoint.hpp"
#include "Core/Event.hpp"

using namespace ComPWA;

double Kinematics::GetPhspVolume() const {
  if (!is_PS_area_calculated_) {
    const_cast<double &>(PS_area_) = calculatePSArea();
    const_cast<bool &>(is_PS_area_calculated_) = true;
  }
  return PS_area_;
}

void Kinematics::SetPhspVolume(double vol) {
  PS_area_ = vol;
  is_PS_area_calculated_ = true;
  
  LOG(info)<<"Kinematics::SetPhspVolume() | Setting phase space "
  "volume to "<<std::to_string(GetPhspVolume())<<".";
}
