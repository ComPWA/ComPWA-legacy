//! Data Interface Base-Class.
/*! \class DIFBase
 * @file DIFBase.hpp
 * This class provides the interface to experimental data. As it is pure virtual,
 * one needs at least one implementation to provide data for the other modules. If
 * a new reader is derived from and fulfills this base-class, no change in other
 * modules are necessary to work with the new dataset.
*/

#ifndef DIFBASE_HPP_
#define DIFBASE_HPP_

#include "PWAEvent.hpp"

class DIFBase
{

public:

  DIFBase()
	  {
	  }

  virtual ~DIFBase()
	{ /* nothing */	}
  virtual const int getEvent(const int, PWAEvent&) =0;

  virtual const unsigned int getNEvents() const =0;

};

#endif
