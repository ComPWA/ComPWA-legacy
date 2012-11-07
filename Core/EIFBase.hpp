//! Estimator Interface Base-Class.
/*! \class EIFBase
 * @file EIFBase.hpp
 * This class provides the interface to classes which estimate the "closeness" of
 * the modeled intensity to the data. As it is pure virtual, one needs at least
 * one implementation to provide an estimator for the analysis. If a new estimator
 * is derived from and fulfills this base-class, no change in other modules are
 * necessary to work with the new estimation function.
*/

#ifndef _EIFBASE_HPP_
#define _EIFBASE_HPP_

class EIFBase
{

public:

  EIFBase(){
  }

  virtual ~EIFBase(){
  /* nothing */
  }

  virtual const int getEstimatedVal(double&) = 0;


};

#endif
