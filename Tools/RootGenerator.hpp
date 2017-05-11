/*
 * RootGenerator.hpp
 *
 *  Created on: 15 Sep 2016
 *      Author: weidenka
 */

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

#ifndef TOOLS_ROOTGENERATOR_HPP_
#define TOOLS_ROOTGENERATOR_HPP_

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <memory>

#include "TLorentzVector.h"
#include "TParticle.h"
#include "TGenPhaseSpace.h"
#include "TRandom3.h"

#include "Core/Generator.hpp"
#include "Core/Event.hpp"
#include "Core/Particle.hpp"
#include "Core/Kinematics.hpp"

namespace ComPWA {
namespace Tools {

class RootGenerator : public Generator {

public:
  //! Constructor for a three particle decay with given masses
  RootGenerator(double sqrtS, double m1, double m2, double m3, int seed = -1);
  
  //! Default Constructor. Information on the decay is obtained from Kinematics
  RootGenerator(int seed = -1);
  
  ~RootGenerator() { delete[] masses; };

  virtual RootGenerator *Clone();
  
  virtual void Generate(Event &evt);
  
  virtual void SetSeed(unsigned int seed);
  
  virtual unsigned int GetSeed() const;
  
  virtual double GetUniform(double min, double max) const;
  
  virtual double GetGaussDist(double mu, double sigma) const;
  
  virtual TGenPhaseSpace *GetGenerator() { return &event; }

protected:
  double sqrtS;
  
  TGenPhaseSpace event;
  
  size_t nPart;
  
  Double_t *masses;
};

class UniformTwoBodyGenerator : public RootGenerator {
public:
  UniformTwoBodyGenerator(double minSq_, double maxSq_, int seed = -1)
      : RootGenerator(seed), minSq(minSq_), maxSq(maxSq_) {
    if (Kinematics::Instance()->GetFinalState().size() != 2)
      throw std::runtime_error("UniformTwoBodyGenerator::"
                               "UniformTwoBodyGenerator() | Not a two body "
                               "decay!");
  }
  virtual void Generate(Event &evt);
  virtual UniformTwoBodyGenerator *Clone() {
    return (new UniformTwoBodyGenerator(*this));
  }

protected:
  double minSq, maxSq;
};

} /* namespace Tools*/
} /* namespace ComPWA */

#endif /* TOOLS_ROOTGENERATOR_HPP_ */
