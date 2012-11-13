//! Optimizer Data Base-Class.
/*! \class OIFData
 * @file OIFData.hpp
 * This class provides the interface for the optimizer module to access the data.
 * If the access to the data is provided fulfilling this interface, then one can
 * use the same data and parameter-set with different optimizers.
*/

#ifndef OIFDATA_HPP_
#define OIFDATA_HPP_

#include <vector>
#include "PWAParameter.hpp"

class OIFData{

public:

  OIFData(){
  }

  virtual ~OIFData(){ /* nothing */	}

  virtual double controlParameter(const std::vector<PWAParameter<double> >& minPar) =0;
 
};

#endif
