//! Optimizer Interface Base-Class.
/*! \class OIFBase
 * @file OIFBase.hpp
 * This class provides the interface to (external) optimization libraries or
 * routines. As it is pure virtual, one needs at least one implementation to
 * provide an optimizer for the analysis which varies free model-parameters. If
 * a new optimizer is derived from and fulfills this base-class, no change in
 * other modules are necessary to work with the new optimizer library or routine.
*/

#ifndef _OIFBASE_HPP_
#define _OIFBASE_HPP_

#include <vector>
#include <memory>

#include "PWAParameter.hpp"

class OIFBase
{

public:

  OIFBase()
	  {
	  }

  virtual ~OIFBase()
	{ /* nothing */	}

  // TODO: template <class T> or map?
  virtual const double exec(std::vector<std::shared_ptr<PWAParameter> >& par) =0;
 
};

#endif
