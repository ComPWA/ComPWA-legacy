// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <algorithm>
#include "Core/Exceptions.hpp"
#include "Core/Kinematics.hpp"
#include "Core/DataPoint.hpp"

namespace ComPWA {

DataPoint::DataPoint() : Weight(1.), Eff(1.) {
  return;
}

void DataPoint::reset(unsigned int size) { Values = std::vector<double>(size); }

void DataPoint::setValue(unsigned int pos, double val) {
  try {
    Values.at(pos) = val;
  } catch (...) {
    LOG(ERROR) << "dataPoint::setVal() | Can not access index " << pos << "!";
    throw;
  }
  return;
}

double DataPoint::value(unsigned int num) const {
  double rt;

  try {
    rt = Values.at(num);
  } catch (...) {
    LOG(ERROR) << "dataPoint::getVal() | Can not access index " << num << "!";
    throw;
  }
  return rt;
}

std::ostream &operator<<(std::ostream &os, const DataPoint &p) {
  os << "(";
  for (int i = 0; i < p.size(); i++){
    os << p.value(i);
  if( i == p.size()-1 )
    os << ")";
  else 
    os << ", ";
  }
  return os;
}

} // ns::ComPWA
