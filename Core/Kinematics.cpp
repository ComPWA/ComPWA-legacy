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

Kinematics* Kinematics::instance() {
  if (!inst_) {
    throw std::runtime_error(
        "No instance of Kinematics created! Create one first!");
  }

  return Kinematics::inst_;
}

Kinematics* Kinematics::inst_ = 0;

Kinematics::Kinematics() {
}

Kinematics::~Kinematics() {
}

