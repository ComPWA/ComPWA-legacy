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

#ifndef GENERATOR_HPP_
#define GENERATOR_HPP_
#include "Core/Event.hpp"

namespace ComPWA {

/**
 *  \class Generator
 *  \brief Virtual class for PHSP generators
 */
class Generator {
private:
public:
  Generator(){};
  virtual ~Generator(){};
  virtual void Generate(Event &) = 0;
  virtual Generator *Clone() = 0;
  virtual void SetSeed(unsigned int) = 0;
  virtual unsigned int GetSeed() const = 0;
  virtual double GetUniform(double min, double max) const = 0;
  virtual double GetGaussDist(double mu, double sigma) const { return 0; }
};

} /* namespace ComPWA */

#endif /* GENERATOR_HPP_ */ 
