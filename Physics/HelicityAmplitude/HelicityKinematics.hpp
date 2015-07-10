//-------------------------------------------------------------------------------
// Copyright (c) 2013 Stefan Pflueger.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//   Stefan Pflueger - initial API and implementation
//-------------------------------------------------------------------------------

#ifndef HELICITYKINEMATICS_HPP_
#define HELICITYKINEMATICS_HPP_

#include "Core/Kinematics.hpp"

namespace HelicityFormalism {

class HelicityKinematics: public Kinematics {
public:
  HelicityKinematics();
  virtual ~HelicityKinematics();

  bool isWithinPhsp(const dataPoint& point) const;
  double getMotherMass() const;
  double getPhspVolume() const;
  void eventToDataPoint(Event& ev, dataPoint& point) const;
  double getMass(unsigned int num) const;
  double getMass(std::string name) const;
};

} /* namespace HelicityFormalism */

#endif /* HELICITYKINEMATICS_HPP_ */
