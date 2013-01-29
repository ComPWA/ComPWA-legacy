//! Data Interface Base-Class.
/*! \class Data
 * @file Data.hpp
 * This class provides the interface to experimental data. As it is pure virtual,
 * one needs at least one implementation to provide data for the other modules. If
 * a new reader is derived from and fulfills this base-class, no change in other
 * modules are necessary to work with the new dataset.
*/

#ifndef DATA_HPP_
#define DATA_HPP_

#include "Core/PWAEvent.hpp"

class Data
{

public:

  Data()
	  {
	  }

  virtual ~Data()
	{ /* nothing */	}
  virtual const int getEvent(const int, PWAEvent&) =0;
  virtual const int getBin(const int, double&, double&)=0; //TODO: seperate?

  virtual const unsigned int getNEvents() const =0;
  virtual const unsigned int getNBins() const =0;

};

#endif
