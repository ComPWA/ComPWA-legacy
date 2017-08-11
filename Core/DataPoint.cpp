// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <algorithm>
#include "Core/Exceptions.hpp"
#include "Core/Kinematics.hpp"
#include "Core/DataPoint.hpp"

namespace ComPWA {

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
