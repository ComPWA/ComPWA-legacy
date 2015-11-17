//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
//
// This file is part of ComPWA, check license.txt for details
//
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//-------------------------------------------------------------------------------
//! Internal container representing a particle.
/*! \class Particle
 * @file Particle.hpp
 * This class provides a internal container for information of a particle. The
 * class provides the momentum 4-vector and pid of the particle.
*/

#ifndef _PARTICLE_HPP_
#define _PARTICLE_HPP_

#include <vector>
#include <math.h>

struct Particle
{

  Particle():px(0),py(0),pz(0),E(0),pid(0){
  }

  Particle(double inPx, double inPy, double inPz, double inE, int inpid=0)
    :px(inPx),py(inPy),pz(inPz),E(inE),pid(inpid){
  }

  //TODO: operators, mass-square of 2 particles ..

  Particle(const Particle& inParticle)
    :px(inParticle.px),py(inParticle.py),pz(inParticle.pz),E(inParticle.E),pid(inParticle.pid){
  }

  inline double invariantMass(const Particle& inParticle) const {
    return (pow(E+inParticle.E,2)-pow(px+inParticle.px ,2)-pow(py+inParticle.py ,2)-pow(pz+inParticle.pz ,2));
  }

  static double invariantMass(const Particle& inPa, const Particle& inPb) {
    return (pow(inPa.E+inPb.E,2)-pow(inPa.px+inPb.px ,2)-pow(inPa.py+inPb.py ,2)-pow(inPa.pz+inPb.pz ,2));
  }

  inline double getMassSquare() const {
    return E*E-px*px-py*py-pz*pz;
  }

 /* virtual const inline double getPx() const {return px;}
  virtual const inline double getPy() const {return py;}
  virtual const inline double getPz() const {return pz;}
  virtual const inline double getE() const {return E;}
  virtual const inline int getPid() const {return pid;}

protected:*/
  double px, py, pz, E;
  int pid;
  //TODO: other particle info?

};

#endif
