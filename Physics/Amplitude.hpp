//! Physics Interface Base-Class.
/*! \class Amplitude
 * @file Amplitude.hpp
 * This class provides the interface to the model which tries to describe the
 * intensity. As it is pure virtual, one needs at least one implementation to
 * provide an model for the analysis which calculates intensities for an event on
 * basis model parameters. If a new physics-model is derived from and fulfills
 * this base-class, no change in other modules are necessary to work with the new
 * physics module.
*/

#ifndef PIFBASE_HPP_
#define PIFBASE_HPP_

#include <vector>
#include <memory>

#include "Core/ParameterList.hpp"

class Amplitude
{

public:

  Amplitude()
	  {
	  }

  virtual ~Amplitude()
	{ /* nothing */	}

  virtual const double integral(ParameterList& par) =0;
  virtual const double volume() =0;

  virtual const double intensity(std::vector<double>& x, ParameterList& par) =0;

  virtual const bool fillStartParVec(ParameterList& outPar) =0;


};

#endif
