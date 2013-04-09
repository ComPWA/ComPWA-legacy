//! Simple \f$\chi^{2}\f$-Estimator.
/*! \class ChiOneD
 * @file ChiOneD.hpp
 * This class calculates a simple \f$\chi^{2}\f$ of a intensity and a dataset.
 * Data and Model are provided in the constructor using the Amplitude and Data
 * interfaces. The class itself fulfills the Estimator interface.
*/

#ifndef _EIFChiOneD_HPP
#define _EIFChiOneD_HPP

#include <vector>
#include <memory>
#include <string>

//PWA-Headers
#include "Estimator/Estimator.hpp"
#include "Physics/Amplitude.hpp"
#include "DataReader/Data.hpp"
#include "Core/Event.hpp"
#include "Core/ParameterList.hpp"

class ChiOneD : public Estimator {

public:
  /// Default Constructor (0x0)
  ChiOneD(std::shared_ptr<Amplitude>, std::shared_ptr<Data>);

  virtual double controlParameter(ParameterList& minPar);

  /** Destructor */
  virtual ~ChiOneD();

protected:

private:
  std::shared_ptr<Amplitude> pPIF_;
  std::shared_ptr<Data> pDIF_;

};

#endif /* _EIFChiOneD_HPP */
