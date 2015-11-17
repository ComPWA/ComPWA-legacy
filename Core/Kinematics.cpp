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

#include <stdexcept>

#include "Core/Kinematics.hpp"
#include "Core/DataPoint.hpp"
#include "Core/Event.hpp"

Kinematics* Kinematics::instance() {
  if (!inst_) {
    throw std::runtime_error(
        "No instance of Kinematics created! Create one first!");
  }

  return Kinematics::inst_;
}

Kinematics* Kinematics::inst_ = 0;

Kinematics::Kinematics() :
    is_PS_area_calculated_(false), PS_area_(0.0) {
}

Kinematics::~Kinematics() {
}

unsigned int Kinematics::getVariableIndex(std::string name) const {
  unsigned int pos(999);
  std::vector<std::string> varNames = Kinematics::instance()->getVarNames();
  std::vector<std::string>::const_iterator search_result = find(varNames.begin(), varNames.end(), name);
  if(search_result != varNames.end()) {
    pos = search_result - varNames.begin();
  }
  else {
    BOOST_LOG_TRIVIAL(error)<<"dataPoint::getVal(): variable with name "<<name<<" not found!";
  }
  return pos;
}

//! calculated the PHSP volume of the current decay by MC integration
double Kinematics::getPhspVolume() {
  if (!is_PS_area_calculated_)
    PS_area_ = calculatePSArea();
  return PS_area_;
}

//! converts Event to dataPoint
void Kinematics::eventToDataPoint(Event& ev, dataPoint& point) const {
  // set event weight as data point weight first
  double weight = ev.getWeight();
  point.setWeight(weight);
  // now do the actual transformation
  translateEventToDataPoint(ev, point);
}
