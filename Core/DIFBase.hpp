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
  virtual const int getEvent(const int, shared_ptr<PWAEvent>) =0;

  virtual const unsigned int getNEvents() const =0;

};

#endif
