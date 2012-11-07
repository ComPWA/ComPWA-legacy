//! Optimizer Interface Base-Class.
/*! \class OIFBase
 * @file OIFBase.hpp
 * This class provides the interface to (external) optimization libraries or
 * routines. As it is pure virtual, one needs at least one implementation to
 * provide an optimizer for the analysis which varies free model-parameters. If
 * a new optimizer is derived from and fulfills this base-class, no change in
 * other modules are necessary to work with the new optimizer library or routine.
*/

#ifndef MIBASE_HPP_
#define MIBASE_HPP_

class OIFBase
{

public:

  OIFBase()
	  {
	  }

  virtual ~OIFBase()
	{ /* nothing */	}

  virtual const double exec(unsigned int Npar, double* par,  double* min, double* max, double* err) =0;
 
};

#endif
