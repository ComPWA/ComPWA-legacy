//-------------------------------------------------------------------------------
// Copyright (c) 2013 Peter Weidenkaff.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Peter Weidenkaff - initial API
//-------------------------------------------------------------------------------

#include <algorithm>
#include "Core/Exceptions.hpp"
#include "Core/Kinematics.hpp"
#include "Core/DataPoint.hpp"

namespace ComPWA {

dataPoint::dataPoint(std::vector<double> vec) : weight(1.), eff(1.) {
  if (Kinematics::Instance()->GetNVars() != vec.size())
    throw std::runtime_error("dataPoint::dataPoint() vector has wrong length!");
  var = vec;
  return;
}

dataPoint::dataPoint(const Event &ev) : weight(1.), eff(1.) {
  Kinematics::Instance()->EventToDataPoint(ev, *this);
  weight = ev.getWeight();
  return;
}

dataPoint::dataPoint() : weight(1.), eff(1.) {
  return;
}

void dataPoint::Reset(unsigned int size) { var = std::vector<double>(size); }

void dataPoint::SetValue(unsigned int pos, double val) {
  try {
    var.at(pos) = val;
  } catch (...) {
    LOG(error) << "dataPoint::setVal() | "
                  "Can not access index "
               << pos << "!";
    throw;
  }
  return;
}

double dataPoint::GetValue(unsigned int num) const {
  double rt;

  try {
    rt = var.at(num);
  } catch (...) {
    LOG(error) << "dataPoint::getVal() | "
                  "Can not access index "
               << num << "!";
    throw;
  }
  return rt;
}

void dataPoint::SetPoint(std::vector<double> values) {
  if (Kinematics::Instance()->GetNVars() != values.size())
    throw std::runtime_error("dataPoint::setPoint() vector has wrong length!");
  var = std::vector<double>(values);
  return;
}

std::ostream &operator<<(std::ostream &os, const dataPoint &p) {
  os << "(";
  for (int i = 0; i < p.Size(); i++){
    os << p.GetValue(i);
  if( i == p.Size()-1 )
    os << ")";
  else 
    os << ", ";
  }
  return os;
}

} /* namespace ComPWA */
