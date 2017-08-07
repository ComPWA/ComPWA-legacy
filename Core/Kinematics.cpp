//-------------------------------------------------------------------------------
// Copyright (c) 2013 Peter Weidenkaff.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//		Peter Weidenkaff -
//-------------------------------------------------------------------------------

#include "Core/Kinematics.hpp"
#include "Core/Exceptions.hpp"
#include "Core/PhysConst.hpp"
#include "Core/DataPoint.hpp"
#include "Core/Event.hpp"

namespace ComPWA {

//! calculated the PHSP volume of the current decay by MC integration
double Kinematics::GetPhspVolume() {
  if (!is_PS_area_calculated_) {
    PS_area_ = calculatePSArea();
    is_PS_area_calculated_ = true;
  }
  return PS_area_;
}

void Kinematics::SetPhspVolume(double vol) {
  PS_area_ = vol;
  is_PS_area_calculated_ = true;
  
  LOG(info)<<"Kinematics::SetPhspVolume() | Setting phase space "
  "volume to "<<std::to_string(GetPhspVolume())<<".";
}

} /* namespace ComPWA */
