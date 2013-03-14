//! Optimizer Data Base-Class.
/*! \class ControlParameter
 * @file ControlParameter.hpp
 * This class provides the interface for the optimizer module to access the data.
 * If the access to the data is provided fulfilling this interface, then one can
 * use the same data and parameter-set with different optimizers.
*/

#ifndef CONTROLPARAMETER_HPP_
#define CONTROLPARAMETER_HPP_

#include <vector>
#include <memory>

#include "Core/ParameterList.hpp"

class ControlParameter{

public:

  ControlParameter(){
  }

  virtual ~ControlParameter(){ /* nothing */	}

  virtual double controlParameter(ParameterList& minPar) =0;
 
};

#endif
