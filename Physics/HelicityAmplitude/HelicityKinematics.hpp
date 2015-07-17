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
  std::vector<std::vector<ParticleState> > decay_trees_;

  bool is_PS_area_calculated_;
  double PS_area_;

  HelicityKinematics(const std::vector<HelicityDecayTree>& decay_trees);
  virtual ~HelicityKinematics();

  // delete methods to ensure that there will only be one instance
  HelicityKinematics(const HelicityKinematics&) = delete;
  void operator=(const HelicityKinematics&) = delete;

  void calculatePSArea();
  Vector4 determineBoostedKinematicVariables(
      std::pair<Vector4, Vector4> two_body_state, Vector4 mother);

public:
  static Kinematics* createInstance(const std::vector<HelicityDecayTree>& decay_trees) {
    if (0 == inst_)
      inst_ = new HelicityKinematics(decay_tree);
    return inst_;
  }

  bool isWithinPhsp(const dataPoint& point) const;
  double getMotherMass() const;
  double getPhspVolume() const;
  void eventToDataPoint(Event& ev, dataPoint& point) const;
  double getMass(unsigned int num) const;
  double getMass(std::string name) const;
};

} /* namespace HelicityFormalism */

#endif /* HELICITYKINEMATICS_HPP_ */
