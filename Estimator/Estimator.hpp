//! Estimator Interface Base-Class.
/*! \class Estimator
 * @file Estimator.hpp
 * This class provides the interface to classes which estimate the "closeness" of
 * the modeled intensity to the data. As it is pure virtual, one needs at least
 * one implementation to provide an estimator for the analysis. If a new estimator
 * is derived from and fulfills this base-class, no change in other modules are
 * necessary to work with the new estimation function. As it is derived from
 * ControlParameter, it can be used in the optimizer modules.
*/

#ifndef _ESTIMATOR_HPP_
#define _ESTIMATOR_HPP_

#include <vector>
#include <memory>

#include "Core/ControlParameter.hpp"
#include "Core/PWAParameter.hpp"

class Estimator : public ControlParameter
{

public:

  Estimator(){
  }

  virtual ~Estimator(){
  /* nothing */
  }

  virtual double controlParameter(std::vector<std::shared_ptr<PWAParameter> >& minPar) = 0;

};

#endif
