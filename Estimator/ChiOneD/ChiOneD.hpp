// Copyright (c) 2013 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

//! Simple \f$\chi^{2}\f$-Estimator.
/*! \class ChiOneD
 * @file ChiOneD.hpp
 * This class calculates a simple \f$\chi^{2}\f$ of a intensity and a dataset.
 * Data and Model are provided in the constructor using the Amplitude and Data
 * interfaces. The class itself fulfills the Estimator interface.
 */

#ifndef COMPWA_ESTIMATOR_CHIONED_HPP_
#define COMPWA_ESTIMATOR_CHIONED_HPP_

#include <memory>
#include <vector>

#include "Estimator/Estimator.hpp"

namespace ComPWA {
namespace Data {
struct DataSet;
}
namespace Estimator {
namespace ChiOneD {

class ChiOneD : public ComPWA::Estimator::Estimator<double> {

public:
  ChiOneD(std::shared_ptr<ComPWA::Intensity> Intensity_,
          const Data::DataSet &DataSample_);

  double evaluate() noexcept;

private:
  std::shared_ptr<ComPWA::Intensity> Intensity;
  const Data::DataSet &DataSample;
};

} // namespace ChiOneD
} // namespace Estimator
} // namespace ComPWA

#endif
