// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Generator base class.
///

#ifndef GENERATOR_HPP_
#define GENERATOR_HPP_

namespace ComPWA {

class Event;

/**
 *  \class Generator
 *  \brief Virtual class for PHSP generators
 */
class Generator {
public:
  virtual void generate(Event &) = 0;
  
  virtual Generator *clone() = 0;
  
  virtual void setSeed(unsigned int) = 0;
  
  virtual unsigned int seed() const = 0;
  
  virtual double uniform(double min, double max) const = 0;
  
  virtual double gauss(double mu, double sigma) const { return 0; }
};

} // ns::ComPWA

#endif
