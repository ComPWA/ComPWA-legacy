// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Generator base class.
///

#ifndef GENERATOR_HPP_
#define GENERATOR_HPP_

#include "Core/Event.hpp"

namespace ComPWA {

/**
 *  \class Generator
 *  \brief Virtual class for PHSP generators
 */
class Generator {
public:
  virtual ~Generator() {};

  virtual ComPWA::Event generate() = 0;
  
  virtual void setSeed(unsigned int) = 0;
  
  virtual unsigned int getSeed() const = 0;
  
  virtual double uniform(double min, double max) = 0;
  
  virtual double gauss(double mu, double sigma) const { return 0; }
};

} // ns::ComPWA

#endif
