//! Negative Log Likelihood-Estimator.
/*! \class EIFMinLogLH
 * @file EIFMinLogLH.hpp
 * This class calculates a simple -log(LH) of a intensity and a dataset.
 * Data and Model are provided in the constructor using the Amplitude and Data
 * interfaces. The class itself fulfills the Estimator interface.
*/

#ifndef _EIFMINLOGLH_HPP
#define _EIFMINLOGLH_HPP

#include <vector>
#include <memory>
#include <string>

//PWA-Headers
#include "Estimator.hpp"
#include "Amplitude.hpp"
#include "Data.hpp"
#include "PWAEvent.hpp"
#include "PWAParameter.hpp"

class EIFMinLogLH : public Estimator {

public:
  /// Default Constructor (0x0)
  EIFMinLogLH(std::shared_ptr<Amplitude>, std::shared_ptr<Data>);

  virtual double controlParameter(std::vector<std::shared_ptr<PWAParameter> >& minPar);

  /** Destructor */
  virtual ~EIFMinLogLH();

protected:

private:
  std::shared_ptr<Amplitude> pPIF_;
  std::shared_ptr<Data> pDIF_;

};

#endif /* _EIFMINLOGLH_HPP */
