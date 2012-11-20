//! Estimator Interface Base-Class.
/*! \class EIFBase
 * @file EIFBase.hpp
 * This class provides the interface to classes which estimate the "closeness" of
 * the modeled intensity to the data. As it is pure virtual, one needs at least
 * one implementation to provide an estimator for the analysis. If a new estimator
 * is derived from and fulfills this base-class, no change in other modules are
 * necessary to work with the new estimation function. As it is derived from
 * OIFData, it can be used in the optimizer modules.
*/

#ifndef _EIFBASE_HPP_
#define _EIFBASE_HPP_

#include "OIFData.hpp"
#include "PWAParameter.hpp"

class EIFBase : public OIFData
{

public:

  EIFBase(){
  }

  virtual ~EIFBase(){
  /* nothing */
  }

  virtual double controlParameter(std::vector<PWAParameter<double> >& minPar) = 0;

};

#endif
