//! Negative Log Likelihood-Estimator.
/*! \class MinLogLH
 * @file MinLogLH.hpp
 * This class calculates a simple -log(LH) of a intensity and a dataset.
 * Data and Model are provided in the constructor using the Amplitude and Data
 * interfaces. The class itself fulfills the Estimator interface.
*/

#ifndef _MINLOGLH_HPP
#define _MINLOGLH_HPP

#include <vector>
#include <memory>
#include <string>

//PWA-Headers
#include "Core/Estimator.hpp"
#include "Core/Amplitude.hpp"
#include "Core/Data.hpp"
#include "Core/PWAEvent.hpp"
#include "Core/PWAParameter.hpp"

class MinLogLH : public Estimator {

public:
  /// Default Constructor (0x0)
  MinLogLH(std::shared_ptr<Amplitude>, std::shared_ptr<Data>);

  virtual double controlParameter(std::vector<std::shared_ptr<PWAParameter> >& minPar);

  /** Destructor */
  virtual ~MinLogLH();

protected:

private:
  std::shared_ptr<Amplitude> pPIF_;
  std::shared_ptr<Data> pDIF_;

};

#endif /* _MINLOGLH_HPP */
