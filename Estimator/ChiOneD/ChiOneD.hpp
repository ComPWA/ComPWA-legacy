//-------------------------------------------------------------------------------
// Copyright (c) 2013 michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     michel - initial API and implementation
//-------------------------------------------------------------------------------
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
#include "Core/Amplitude.hpp"
#include "DataReader/Data.hpp"
#include "Core/Event.hpp"
#include "Core/ParameterList.hpp"

namespace ComPWA {
namespace Estimator {
namespace ChiOneD {

class ChiOneD : public Estimator {

public:
  static std::shared_ptr<ComPWA::ControlParameter> createInstance(std::shared_ptr<Amplitude>, std::shared_ptr<DataReader::Data>);
  virtual double controlParameter(ParameterList& minPar);

  /** Destructor */
  virtual ~ChiOneD();

protected:
  /// Default Constructor (0x0)
  ChiOneD(std::shared_ptr<Amplitude>, std::shared_ptr<DataReader::Data>);

private:
  std::shared_ptr<Amplitude> pPIF_;
  std::shared_ptr<DataReader::Data> pDIF_;

};

} /* namespace ChiOneD */
} /* namespace Estimator */
} /* namespace ComPWA */

#endif /* _EIFChiOneD_HPP */
