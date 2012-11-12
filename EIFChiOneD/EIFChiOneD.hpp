//! Simple \f$\chi^{2}\f$-Estimator.
/*! \class EIFChiOneD
 * @file EIFChiOneD.hpp
 * This class calculates a simple \f$\chi^{2}\f$ of a intensity and a dataset.
 * Data and Model are provided in the constructor using the PIFBase and DIFBase
 * interfaces. The class itself fulfills the EIFBase interface.
*/

#ifndef _EIFChiOneD_HPP
#define _EIFChiOneD_HPP

#include <vector>
#include <memory>
#include <string>

//PWA-Headers
#include "EIFBase.hpp"
#include "PIFBase.hpp"
#include "DIFBase.hpp"
#include "PWAEvent.hpp"

class EIFChiOneD : public EIFBase {

public:
  /// Default Constructor (0x0)
  EIFChiOneD(std::shared_ptr<PIFBase>, std::shared_ptr<DIFBase>);

  virtual double controlParameter(const std::vector<double>& minPar);

  /** Destructor */
  virtual ~EIFChiOneD();

protected:

private:
  std::shared_ptr<PIFBase> pPIF_;
  std::shared_ptr<DIFBase> pDIF_;

};

#endif /* _EIFChiOneD_HPP */
