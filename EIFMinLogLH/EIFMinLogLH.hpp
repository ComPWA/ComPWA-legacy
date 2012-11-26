//! Negative Log Likelihood-Estimator.
/*! \class EIFMinLogLH
 * @file EIFMinLogLH.hpp
 * This class calculates a simple -log(LH) of a intensity and a dataset.
 * Data and Model are provided in the constructor using the PIFBase and DIFBase
 * interfaces. The class itself fulfills the EIFBase interface.
*/

#ifndef _EIFMINLOGLH_HPP
#define _EIFMINLOGLH_HPP

#include <vector>
#include <memory>
#include <string>

//PWA-Headers
#include "EIFBase.hpp"
#include "PIFBase.hpp"
#include "DIFBase.hpp"
#include "PWAEvent.hpp"
#include "PWAParameter.hpp"

class EIFMinLogLH : public EIFBase {

public:
  /// Default Constructor (0x0)
  EIFMinLogLH(std::shared_ptr<PIFBase>, std::shared_ptr<DIFBase>);

  virtual double controlParameter(std::vector<std::shared_ptr<PWAParameter> >& minPar);

  /** Destructor */
  virtual ~EIFMinLogLH();

protected:

private:
  std::shared_ptr<PIFBase> pPIF_;
  std::shared_ptr<DIFBase> pDIF_;

};

#endif /* _EIFMINLOGLH_HPP */
