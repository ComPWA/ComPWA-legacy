//-------------------------------------------------------------------------------
// Copyright (c) 2013 Peter Weidenkaff.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//    Peter Weidenkaff -
//-------------------------------------------------------------------------------

#ifndef TWOBODYKINEMATICS_HPP_
#define TWOBODYKINEMATICS_HPP_

#include "Core/Kinematics.hpp"

class TwoBodyKinematics: public Kinematics {
public:
  TwoBodyKinematics(std::string _nameMother, std::string _name1,
      std::string _name2, double deltaMassWindow = 0.5);
  void init();
  static Kinematics* createInstance(std::string _nameMother, std::string _name1,
      std::string _name2, double massWindow = 0.5) {
    if (0 == inst_)
      inst_ = new TwoBodyKinematics(_nameMother, _name1, _name2, massWindow);
    return inst_;
  }
  //! checks of data point is within phase space boundaries
  virtual bool isWithinPhsp(const dataPoint& point);
  //! mass of mother particle
  virtual double getMotherMass() const {
    return M;
  }
  //! calculate phase space area with simple interval
  virtual double calculatePSArea() {
    return (mass_max - mass_min);
  }
  //! converts Event to dataPoint
  void translateEventToDataPoint(const Event& ev, dataPoint& point) const;
  //! get mass of particles
  virtual double getMass(unsigned int num) const;
  //! get mass of paticles
  virtual double getMass(std::string name) const;

protected:
  std::string nameMother;    //! name of mother particle
  double Msq;    //! mass squared of mother particle
  double M;    //! mass of mother particle
  unsigned int spinM;    //! spin of mother particle
  double Br;    //! width of decaying particle

  std::string name1;    //! name of daughter 1
  double mSq1;    //! masse squared of daughter 1
  double m1;    //! masses of daughter 1
  unsigned int spin1;    //! spin of daughter 1
  std::string name2;    //! name of daughter 2
  double mSq2;    //! masse squared of daughter 2
  double m2;    //! masses of daughter 2
  unsigned int spin2;    //! spin of daughter 2

  double mass_sq_min;    //!minimum value of masssq
  double mass_sq_max;    //!maximum value of masssq
  double mass_min;    //!minimum value of masssq
  double mass_max;    //!maximum value of masssq

};

#endif /* TWOBODYKINEMATICS_HPP_ */
